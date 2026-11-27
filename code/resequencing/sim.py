"""Stochastic route execution and SAA resequencing oracle for Gate 1.

Conventions
-----------
A route is given in *sequence order* (index 0 = station, then stops in the
dispatched/actual order). At a check-in after m served stops, the vehicle is at
position m. Local indices for the remaining problem: 0 = current position,
1..U = remaining stops in their current order.

All times are seconds relative to the route departure.
"""
from __future__ import annotations

import numpy as np
from numba import njit

INF = 1e18


# ----------------------------------------------------------------------------
# Scenario generation (numpy, outside numba)
# ----------------------------------------------------------------------------

def posterior_theta(obs_log_ratio, sR, sE):
    """Normal-normal posterior of the route shock given observed leg log-ratios.

    obs_log_ratio = log(T/tbar) + (sR^2+sE^2)/2 for completed legs = theta + eps.
    """
    m = len(obs_log_ratio)
    if m == 0:
        return 0.0, sR ** 2
    v = 1.0 / (1.0 / sR ** 2 + m / sE ** 2)
    mu = v * np.sum(obs_log_ratio) / sE ** 2
    return mu, v


def draw_future(rng, tbar_loc, sbar_loc, theta_mu, theta_var, sR, sE, sS, n_scen,
                phi_mu=0.0, phi_var=0.0, sP=0.0):
    """Draw travel-time matrices and service times for the remaining problem.

    tbar_loc: (N, N) mean travel times among local nodes (0 = current position).
    sbar_loc: (N,) mean service times (entry 0 unused).
    theta: route-level travel shock (posterior N(theta_mu, theta_var)).
    phi: route-level service shock (posterior N(phi_mu, phi_var)); sP = 0 disables it
         (pre-registered setting).
    Returns T (S, N, N) and Sv (S, N).
    """
    N = tbar_loc.shape[0]
    theta = rng.normal(theta_mu, np.sqrt(theta_var), size=n_scen)
    eps = rng.normal(0.0, sE, size=(n_scen, N, N))
    corr = (sR ** 2 + sE ** 2) / 2.0
    T = tbar_loc[None, :, :] * np.exp(theta[:, None, None] + eps - corr)
    zeta = rng.normal(0.0, sS, size=(n_scen, N))
    if sP > 0:
        phi = rng.normal(phi_mu, np.sqrt(phi_var), size=n_scen)
        Sv = sbar_loc[None, :] * np.exp(phi[:, None] + zeta - (sP ** 2 + sS ** 2) / 2.0)
    else:
        Sv = sbar_loc[None, :] * np.exp(zeta - sS ** 2 / 2.0)
    Sv[:, 0] = 0.0
    return T, Sv


def realize_route(rng, tbar_seq, sbar_seq, sR, sE, sS, sP=0.0):
    """Realize one full execution of the dispatched sequence.

    Returns realized leg times (n,), service times (n,), theta and phi.
    Leg k goes from position k-1 to k (leg 0 unused).
    """
    n = len(sbar_seq)
    theta = rng.normal(0.0, sR)
    corr = (sR ** 2 + sE ** 2) / 2.0
    legs = np.zeros(n)
    for k in range(1, n):
        legs[k] = tbar_seq[k - 1, k] * np.exp(theta + rng.normal(0.0, sE) - corr)
    zeta = rng.normal(0.0, sS, size=n)
    if sP > 0:
        phi = rng.normal(0.0, sP)
        svc = sbar_seq * np.exp(phi + zeta - (sP ** 2 + sS ** 2) / 2.0)
    else:
        phi = 0.0
        svc = sbar_seq * np.exp(zeta - sS ** 2 / 2.0)
    svc[0] = 0.0
    return legs, svc, theta, phi


# ----------------------------------------------------------------------------
# Lateness evaluation (numba)
# ----------------------------------------------------------------------------

@njit(cache=True)
def eval_seq(seq, t0, T, Sv, opn, cls, w, wait, cnt=0.0):
    """Mean total weighted lateness (seconds) of visiting `seq` (local ids 1..U);
    with cnt = 1.0, the mean weighted number of late stops instead.

    seq: int64 array of local ids; start at node 0 at time t0.
    T: (S, N, N); Sv: (S, N); opn/cls/w: (N,) per local node.
    wait: if True, arrival before opn waits until opn (TW regime).
    """
    S = T.shape[0]
    total = 0.0
    for s in range(S):
        t = t0
        prev = 0
        late = 0.0
        for k in range(seq.shape[0]):
            j = seq[k]
            t += T[s, prev, j]
            if wait and t < opn[j]:
                t = opn[j]
            d = t - cls[j]
            if d > 0.0:
                late += w[j] * (cnt + (1.0 - cnt) * d)
            t += Sv[s, j]
            prev = j
        total += late
    return total / S


@njit(cache=True)
def _prefix_states(seq, t0, T, Sv, opn, cls, w, wait, pt, pl, cnt):
    """Fill prefix arrays: pt[k, s] time after serving position k-1 (pt[0]=t0),
    pl[k, s] accumulated lateness after k stops."""
    S = T.shape[0]
    U = seq.shape[0]
    for s in range(S):
        pt[0, s] = t0
        pl[0, s] = 0.0
        t = t0
        prev = 0
        late = 0.0
        for k in range(U):
            j = seq[k]
            t += T[s, prev, j]
            if wait and t < opn[j]:
                t = opn[j]
            d = t - cls[j]
            if d > 0.0:
                late += w[j] * (cnt + (1.0 - cnt) * d)
            t += Sv[s, j]
            prev = j
            pt[k + 1, s] = t
            pl[k + 1, s] = late


@njit(cache=True)
def _eval_from(cand, a, seq_prev_node, T, Sv, opn, cls, w, wait, pt, pl, bound, cnt):
    """Evaluate candidate sequence reusing prefix up to position a (exclusive).

    Returns total lateness summed over scenarios, or a value >= bound if pruned.
    seq_prev_node: local id of the node at position a-1 (0 if a == 0).
    """
    S = T.shape[0]
    U = cand.shape[0]
    # running per-scenario state
    tot = 0.0
    for s in range(S):
        tot += pl[a, s]
    if tot >= bound:
        return tot
    t = np.empty(S)
    late = np.empty(S)
    for s in range(S):
        t[s] = pt[a, s]
        late[s] = pl[a, s]
    prev = seq_prev_node
    for k in range(a, U):
        j = cand[k]
        tot = 0.0
        for s in range(S):
            ts = t[s] + T[s, prev, j]
            if wait and ts < opn[j]:
                ts = opn[j]
            d = ts - cls[j]
            if d > 0.0:
                late[s] += w[j] * (cnt + (1.0 - cnt) * d)
            t[s] = ts + Sv[s, j]
            tot += late[s]
        prev = j
        if tot >= bound:
            return tot
    return tot


@njit(cache=True)
def local_search(seq0, t0, T, Sv, opn, cls, w, wait, max_passes, cnt=0.0):
    """First-improvement Or-opt (segments 1-3) + 2-opt descent on the SAA objective.

    Returns (best_seq, best_mean_cost, passes).
    """
    S = T.shape[0]
    U = seq0.shape[0]
    seq = seq0.copy()
    pt = np.empty((U + 1, S))
    pl = np.empty((U + 1, S))
    _prefix_states(seq, t0, T, Sv, opn, cls, w, wait, pt, pl, cnt)
    best = 0.0
    for s in range(S):
        best += pl[U, s]
    if best <= 0.0 or U < 2:
        return seq, best / S, 0
    cand = np.empty(U, dtype=np.int64)
    eps = 1e-9 * S if cnt == 0.0 else 0.5
    passes = 0
    improved = True
    while improved and passes < max_passes:
        improved = False
        passes += 1
        # ---- Or-opt: move segment [i, i+L) to position before index p (in reduced list)
        for L in range(1, 4):
            if improved:
                break
            for i in range(0, U - L + 1):
                if improved:
                    break
                for p in range(0, U - L + 1):
                    if p == i:
                        continue
                    # build candidate: remove segment, insert at p of remaining list
                    c = 0
                    r = 0
                    for q in range(U):
                        if q >= i and q < i + L:
                            continue
                        if r == p:
                            for z in range(L):
                                cand[c] = seq[i + z]
                                c += 1
                        cand[c] = seq[q]
                        c += 1
                        r += 1
                    if r == p:
                        for z in range(L):
                            cand[c] = seq[i + z]
                            c += 1
                    a = i if i < p else p
                    prev_node = 0 if a == 0 else seq[a - 1]
                    val = _eval_from(cand, a, prev_node, T, Sv, opn, cls, w, wait, pt, pl, best - eps, cnt)
                    if val < best - eps:
                        for q in range(U):
                            seq[q] = cand[q]
                        best = val
                        _prefix_states(seq, t0, T, Sv, opn, cls, w, wait, pt, pl, cnt)
                        improved = True
                        break
        if improved:
            continue
        # ---- 2-opt: reverse seq[i..j]
        for i in range(0, U - 1):
            if improved:
                break
            for j in range(i + 1, U):
                for q in range(U):
                    cand[q] = seq[q]
                for q in range(0, j - i + 1):
                    cand[i + q] = seq[j - q]
                prev_node = 0 if i == 0 else seq[i - 1]
                val = _eval_from(cand, i, prev_node, T, Sv, opn, cls, w, wait, pt, pl, best - eps, cnt)
                if val < best - eps:
                    for q in range(U):
                        seq[q] = cand[q]
                    best = val
                    _prefix_states(seq, t0, T, Sv, opn, cls, w, wait, pt, pl, cnt)
                    improved = True
                    break
    return seq, best / S, passes
