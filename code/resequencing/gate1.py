"""Gate 1: recoverable vs. locked-in lateness on Amazon routes (see analysis/gate1/preregistration.md).

Usage:
    python -I code/resequencing/gate1.py <routes.pkl> <out.csv> [--noise base,low,high]
           [--fracs 0.35,0.5] [--regimes ETA,TW] [--sample N] [--workers W] [--seed 2026]
"""
from __future__ import annotations

import argparse
import os
import pickle
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim import draw_future, eval_seq, local_search, posterior_theta, realize_route  # noqa: E402

NOISE = {"low": (0.10, 0.15, 0.25), "base": (0.15, 0.25, 0.40), "high": (0.25, 0.35, 0.60)}
ETA_HALF_WIDTH = 3600.0
S_IN, S_OOS, MAX_PASSES = 48, 256, 60

_ROUTES = None


def _init(path):
    global _ROUTES
    with open(path, "rb") as fh:
        _ROUTES = pickle.load(fh)


def nominal_eta(r):
    """Planned arrival times along the dispatched (actual) sequence with mean times, no waiting."""
    n = r["n"]
    eta = np.zeros(n)
    t = 0.0
    for k in range(1, n):
        t += float(r["tbar"][k - 1, k])
        eta[k] = t
        t += r["sbar"][k]
    return eta


def nominal_schedule_wait(r, tw_open=None):
    """Planned service-start times along the dispatched sequence with mean times, waiting for
    customer windows to open. Returns (start (n,), planned completion of the last stop)."""
    n = r["n"]
    opn = r["tw_open"] if tw_open is None else tw_open
    start = np.zeros(n)
    t = 0.0
    for k in range(1, n):
        t += float(r["tbar"][k - 1, k])
        o = opn[k]
        if not np.isnan(o) and t < o:
            t = o
        start[k] = t
        t += r["sbar"][k]
    return start, t


# ---- Gate 2c dose regimes: day-definite promises plus a target share of customer-window stops
DOSE_SHARE = {"DAY": 0.0, "MIXB20": 0.20, "MIXB50": 0.50}
_OFFSETS = None


def _offsets():
    """Empirical (open, close) offsets of observed windows relative to the planned no-wait ETA,
    fixed from the evaluation split before Gate 2c (analysis/gate2c/window_offsets_eval.npz)."""
    global _OFFSETS
    if _OFFSETS is None:
        here = os.path.dirname(os.path.abspath(__file__))
        z = np.load(os.path.join(here, "..", "..", "analysis", "gate2c", "window_offsets_eval.npz"))
        _OFFSETS = (z["open_off"], z["close_off"])
    return _OFFSETS


def regime_windows(r, regime, eta, seed, ridx):
    """Customer-window arrays (full route) for a regime. Observed windows for ETA/TW/MIXA/MIXB;
    none for DAY; observed + synthetic windows up to the target share for MIXB20/MIXB50.
    The synthetic draw depends only on (seed, route, regime), so it is identical across noise
    levels and check-ins (common random numbers)."""
    o, c = r["tw_open"].copy(), r["tw_close"].copy()
    if regime == "DAY":
        o[:] = np.nan
        c[:] = np.nan
    elif regime in ("MIXB20", "MIXB50"):
        share = DOSE_SHARE[regime]
        n = r["n"]
        free = np.where(np.isnan(c[1:]))[0] + 1
        need = max(0, int(round(share * (n - 1))) - int(np.sum(~np.isnan(c[1:]))))
        if need > 0 and len(free) > 0:
            rng = np.random.default_rng([seed, ridx, 9, int(share * 100)])
            pos = rng.choice(free, size=min(need, len(free)), replace=False)
            op, cl = _offsets()
            j = rng.integers(0, len(op), size=len(pos))
            o[pos] = np.where(op[j] <= -1e8, np.nan, eta[pos] + op[j])
            c[pos] = eta[pos] + cl[j]
    return o, c


def promise_arrays(r, idx, regime, eta, eta_hw, day_deadline, start_w=None, tw_open=None, tw_close=None):
    """Opening times, deadlines, weights and waiting flag for local nodes idx under a promise regime.

    MIXA uses the planned schedule *with waiting* (start_w) as the ETA, so no stop is promised
    before its window can open. DAY / MIXB20 / MIXB50 are MIXB with modified window arrays."""
    if regime in DOSE_SHARE:
        regime = "MIXB"
    tw_o = (r["tw_open"] if tw_open is None else tw_open)[idx]
    tw_c = (r["tw_close"] if tw_close is None else tw_close)[idx]
    has_tw = ~np.isnan(tw_c)
    if regime == "TW":
        opn = np.where(np.isnan(tw_o), -1e18, tw_o)
        cls = np.where(has_tw, tw_c, 1e18)
        w = has_tw.astype(np.float64)
        wait = True
    elif regime == "MIXB":  # day-definite + customer windows
        opn = np.where(np.isnan(tw_o), -1e18, tw_o)
        cls = np.where(has_tw, tw_c, day_deadline)
        w = np.ones(len(idx))
        wait = True
    elif regime == "MIXA":  # communicated ETA + customer windows
        opn = np.where(np.isnan(tw_o), -1e18, tw_o)
        base = start_w[idx] + eta_hw
        cls = np.where(has_tw, np.minimum(np.where(has_tw, tw_c, 1e18), base), base)
        w = np.ones(len(idx))
        wait = True
    else:  # ETA
        opn = np.full(len(idx), -1e18)
        cls = eta[idx] + eta_hw
        w = np.ones(len(idx))
        wait = False
    w = w.copy()
    w[0] = 0.0
    return opn, cls, w, wait


def run_route(task):
    ridx, noise, fracs, regimes, seed, eta_hw, sP, cnt = task
    r = _ROUTES[ridx]
    sR, sE, sS = NOISE[noise]
    n = r["n"]
    tbar = r["tbar"].astype(np.float64)
    sbar = r["sbar"]
    rng_real = np.random.default_rng([seed, ridx, 0, list(NOISE).index(noise)])
    legs, svc, theta_true, phi_true = realize_route(rng_real, tbar, sbar, sR, sE, sS, sP)
    eta = nominal_eta(r)
    # per-regime window arrays and planned schedules (with waiting) -> day-definite deadline
    WIN = {}
    for regime in regimes:
        o_full, c_full = regime_windows(r, regime, eta, seed, ridx)
        st_w, en_w = nominal_schedule_wait(r, o_full)
        WIN[regime] = (o_full, c_full, st_w, en_w + eta_hw)
    corr = (sR ** 2 + sE ** 2) / 2.0
    corr_s = (sP ** 2 + sS ** 2) / 2.0
    out = []
    for f in fracs:
        m = int(round(f * (n - 1)))
        if m < 1 or n - 1 - m < 2:
            continue
        idx = np.arange(m, n)  # local 0 = current position m
        U = len(idx) - 1
        tbar_loc = tbar[np.ix_(idx, idx)]
        sbar_loc = sbar[idx].copy()
        obs = [np.log(legs[k] / tbar[k - 1, k]) + corr for k in range(1, m + 1) if tbar[k - 1, k] > 0]
        mu, var = posterior_theta(np.array(obs), sR, sE)
        if sP > 0:
            obs_s = [np.log(svc[k] / sbar[k]) + corr_s for k in range(1, m + 1) if sbar[k] > 0]
            pmu, pvar = posterior_theta(np.array(obs_s), sP, sS)
        else:
            pmu, pvar = 0.0, 0.0
        rng_in = np.random.default_rng([seed, ridx, 1, list(NOISE).index(noise), int(f * 100)])
        rng_oos = np.random.default_rng([seed, ridx, 2, list(NOISE).index(noise), int(f * 100)])
        T_in, Sv_in = draw_future(rng_in, tbar_loc, sbar_loc, mu, var, sR, sE, sS, S_IN, pmu, pvar, sP)
        T_oos, Sv_oos = draw_future(rng_oos, tbar_loc, sbar_loc, mu, var, sR, sE, sS, S_OOS, pmu, pvar, sP)
        # Gate 2 information ladder: prior-only forecast (route shock not learned) on its own stream
        rng_prior = np.random.default_rng([seed, ridx, 3, list(NOISE).index(noise), int(f * 100)])
        T_pr, Sv_pr = draw_future(rng_prior, tbar_loc, sbar_loc, 0.0, sR ** 2, sR, sE, sS, S_OOS, 0.0, sP ** 2, sP)
        T_det, Sv_det = tbar_loc[None, :, :].copy(), sbar_loc[None, :].copy()
        Sv_det[:, 0] = 0.0
        for regime in regimes:
            o_full, c_full, start_w, day_deadline = WIN[regime]
            opn, cls, w, wait = promise_arrays(r, idx, regime, eta, eta_hw, day_deadline, start_w, o_full, c_full)
            # time at check-in: realized execution of the first m stops under this regime's waiting rule
            t0 = 0.0
            for k in range(1, m + 1):
                t0 += legs[k]
                if wait and not np.isnan(o_full[k]) and t0 < o_full[k]:
                    t0 = o_full[k]
                t0 += svc[k]
            cur = np.arange(1, U + 1, dtype=np.int64)
            tic = time.perf_counter()
            rho_in = eval_seq(cur, t0, T_in, Sv_in, opn, cls, w, wait, cnt)
            best, best_in, passes = cur, rho_in, 0
            if rho_in > 0:
                s1, c1, p1 = local_search(cur, t0, T_in, Sv_in, opn, cls, w, wait, MAX_PASSES, cnt)
                best, best_in, passes = s1, c1, p1
                if regime != "ETA":
                    edd = (np.argsort(cls[1:], kind="stable") + 1).astype(np.int64)
                    if not np.array_equal(edd, cur):
                        s2, c2, p2 = local_search(edd, t0, T_in, Sv_in, opn, cls, w, wait, MAX_PASSES, cnt)
                        if c2 < best_in:
                            best, best_in, passes = s2, c2, passes + p2
            rho = eval_seq(cur, t0, T_oos, Sv_oos, opn, cls, w, wait, cnt)
            lam = eval_seq(best, t0, T_oos, Sv_oos, opn, cls, w, wait, cnt) if passes > 0 else rho
            moved = int(np.sum(best != cur))
            scale = 1.0 if cnt else 60.0  # count objective: values are expected numbers of late stops
            # Gate 2 ladder scores for the current sequence (all on the same promise arrays)
            L2 = eval_seq(cur, t0, T_det, Sv_det, opn, cls, w, wait, cnt)
            L3 = eval_seq(cur, t0, T_pr, Sv_pr, opn, cls, w, wait, cnt)
            b_cnt = eval_seq(cur, t0, T_oos, Sv_oos, opn, cls, w, wait, 1.0)
            has_win = ~np.isnan(c_full[idx])
            has_win[0] = False
            if regime == "ETA":  # windows are not promises under ETA promises
                has_win[:] = False
            b_win = eval_seq(cur, t0, T_oos, Sv_oos, opn, cls, w * has_win, wait, 1.0) if has_win.any() else 0.0
            plan_done = (start_w[m] if wait else eta[m]) + sbar[m]
            out.append(dict(
                route=r["id"], ridx=ridx, station=r["station"], n=n, U=U, frac=f, regime=regime, noise=noise,
                eta_hw_min=eta_hw / 60, sP=sP, objective=('count' if cnt else 'minutes'),
                t0_min=t0 / 60, delay_so_far_min=(t0 - (eta[m] + sbar[m])) / 60, theta_true=theta_true, theta_post=mu,
                phi_true=phi_true, phi_post=pmu,
                tw_remaining=int(np.sum(w[1:] > 0)),
                rho_min=rho / scale, lam_min=lam / scale, delta_min=(rho - lam) / scale,
                rho_in_min=rho_in / scale, best_in_min=best_in / scale, passes=passes, stops_moved=moved,
                delay_plan_min=(t0 - plan_done) / 60, L2_proj=L2 / scale, L3_prior=L3 / scale, b_late_stops=b_cnt,
                b_late_window_stops=b_win, windows_remaining=int(np.sum(has_win)),
                window_share_route=float(np.mean(~np.isnan(c_full[1:]))),
                secs=time.perf_counter() - tic,
            ))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("routes")
    ap.add_argument("out")
    ap.add_argument("--noise", default="base")
    ap.add_argument("--fracs", default="0.35,0.5")
    ap.add_argument("--regimes", default="ETA,TW")
    ap.add_argument("--sample", type=int, default=0)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--eta-half-width-min", type=float, default=ETA_HALF_WIDTH / 60,
                    help="ETA promise = ETA + this many minutes (pre-registered: 60)")
    ap.add_argument("--svc-route-shock", type=float, default=0.0,
                    help="sd of route-level service-time shock (pre-registered: 0 = off)")
    ap.add_argument("--objective", default="minutes", choices=["minutes", "count"],
                    help="minutes = total lateness (pre-registered); count = number of late stops")
    a = ap.parse_args()
    with open(a.routes, "rb") as fh:
        routes = pickle.load(fh)
    unknown = set(a.regimes.split(",")) - {"ETA", "TW", "MIXA", "MIXB"} - set(DOSE_SHARE)
    if unknown:
        raise SystemExit(f"unknown regimes: {unknown}")
    elig = [i for i, r in enumerate(routes) if r["n"] - 1 >= 10]
    if a.sample and a.sample < len(elig):
        rng = np.random.default_rng(a.seed)
        stations = pd.Series([routes[i]["station"] for i in elig])
        frac = a.sample / len(elig)
        pick = []
        for st, grp in stations.groupby(stations):
            k = max(1, int(round(frac * len(grp))))
            pick += list(rng.choice(np.array(elig)[grp.index], size=min(k, len(grp)), replace=False))
        elig = sorted(int(i) for i in pick)
    del routes
    fracs = [float(x) for x in a.fracs.split(",")]
    regimes = a.regimes.split(",")
    tasks = [(i, nz, fracs, regimes, a.seed, a.eta_half_width_min * 60.0, a.svc_route_shock,
              1.0 if a.objective == "count" else 0.0)
             for nz in a.noise.split(",") for i in elig]
    print(f"routes {len(elig)}; tasks {len(tasks)}; workers {a.workers}", flush=True)
    rows = []
    t0 = time.time()
    with Pool(a.workers, initializer=_init, initargs=(a.routes,)) as pool:
        for k, res in enumerate(pool.imap_unordered(run_route, tasks, chunksize=2)):
            rows += res
            if (k + 1) % 100 == 0:
                print(f"{k + 1}/{len(tasks)} tasks, {time.time() - t0:.0f}s", flush=True)
                pd.DataFrame(rows).to_csv(a.out, index=False)
    pd.DataFrame(rows).to_csv(a.out, index=False)
    print(f"done in {time.time() - t0:.0f}s -> {a.out}")


if __name__ == "__main__":
    main()
