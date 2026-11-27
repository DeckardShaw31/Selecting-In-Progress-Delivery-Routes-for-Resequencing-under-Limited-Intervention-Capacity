"""Numerical checks of the theory results, numbered as in manuscript.tex (Lemmas 1-7, Propositions 1-7,
Corollaries 1-4, Remark 1). analysis/theory/proof_audit.md sections 1-3c used an older numbering; the map is in
replication/README.md.

T1  Representation: g_r(x) = integral_{-inf}^x G_r(y) dy, G_r = weighted CDF of effective slack Y = s - xi.
T2  Proposition 2 and Corollary 1: loss <= 2 min(K, N-K) eps*, eps* = 1/2 max_{D_a <= D_b} (rho_a - rho_b)^+,
    and the Wasserstein form loss <= 2 min(K, N-K) W max_r W1(nu_r, nu_bar) (equal total weights W).
T3  Proposition 3, worst case at K = 1: sup over delay vectors of the loss equals Gamma = sup_x [max_r g_r - min_r g_r]
    (brute force over a delay grid must not exceed Gamma and must come close to it).
T4  Proposition 5 (route-level step): |rho_r - pi_r| <= tau_r for any deterministic projection.
T5  Proposition 4 (tail-average bound): loss <= T_m(e+) + T_m(e-) <= sqrt(2m) ||e||_2 for any nondecreasing gbar
    (isotonic and random monotone fits); tail bound <= the Proposition 2 bound for the same gbar;
    LP-optimized tail bound between the loss and the isotonic tail bound; isotonic fit minimizes ||e||_2.
T6  Proposition 4 tightness: the construction in the proof attains both inequalities (K <= N-K and K > N-K).
T7  Corollary 2: |g_r - gbar| <= W d_r with d_r = W1(nu_r, nu_bar), and loss <= W T_2m(d) <= W sqrt(2m) ||d||_2.
T8  Lemma 3 (linear learning): in the linear-Gaussian route-shock model, Monte Carlo risk at a given delay equals
    the closed form with lambda = 1 + c beta and Gaussian residual, and equals the integral of G-tilde.
T9  Proposition 6 (signal-to-heterogeneity): L / V(S_rho) <= T_2m(d~) / (K E_K) <= 2 (m/K) max d~ / E_K under
    Assumption 2 with equal delay sensitivity.
T10 Proposition 7: beta increasing in phi; Var D increasing with travel-part ratio n_eff (phi = 1 vs 0);
    exact lognormal Var D formula against Monte Carlo; Var(theta | D) formula against Monte Carlo.

Added 2026-11-26 with the proof expansion (manuscript numbering):
T11 Lemma 2 (tail top-K) for ARBITRARY scores: loss <= T_m(e+) + T_m(e-) <= T_2m(|e|) <= 2m max|e|, and the
    tail form of Proposition 5: loss of the projection <= T_2m(tau) for any tau_r >= |rho_r - pi_r|.
T12 Proposition 2 attained for every N and K (three-group instances in the appendix proof).
T13 Propositions 1 and 3 on a BOUNDED delay interval X, any K: the appendix constructions (other routes tied at
    the lower delay, no tie-breaking needed) give S_D = T + {r} and a loss within 1e-6 of the target gap.
T14 Corollary 3 (large fleets): (x - t)^+ <= x^2 / (4t); Tbar_alpha(X) <= sqrt(alpha E X^2); and on a large
    i.i.d. fleet the loss per route is below Tbar_alpha(e+) + Tbar_alpha(e-) for a fixed nondecreasing g.
T15 Lemma 5 (size of prediction errors): E|X - 1| = 2 (2 Phi(s/2) - 1) for X = exp(s Z - s^2/2), Monte Carlo.
T16 Lemma 6 (random effects): E MSW = psi_leg^2 and E MSB = psi_leg^2 + J0 psi_route^2 on unbalanced groups,
    Monte Carlo; and the estimator equals gate3_calibrate.re_decomp.

Added 2026-11-27:
T17 Proposition 5 tight: two-group instances with deterministic arrivals and no ties in pi give
    T_2m(tau) = 2m max tau and a loss of m (2c - eps), i.e. a share 1 - eps / (2c) of the bound.
T18 Corollary 4 (projection from conditional means): on random fleets whose arrival times follow recursion (9)
    with waiting, under the empirical law of M scenarios (so every expectation is exact), the schedule from the
    mean times satisfies Ahat <= E A, hence pi <= rho, and the loss of the projection is <= T_m(tau) for every K;
    the tightness construction in the appendix gives a loss share 1 - eps / c of T_m(tau) = m max tau.
(T10 already checks Lemma 7, the exact lognormal Var D formula, against Monte Carlo.)

Added 2026-11-27 (fourth review pass; Lemma 4 is new, so the appendix-B lemmas are now Lemmas 5-7):
T19 Lemma 4 (dispersion term under Gaussian deviations): with A - Ahat ~ N(0, sigma^2) and projected slack z,
    the term of tau_r is sigma phi(z) for z >= 0 and sigma sqrt(2/pi) for z < 0, and the error of the projected
    lateness is sigma L(|z|), L(u) = phi(u) - u (1 - Phi(u)); Monte Carlo, plus 0 <= L <= phi and L decreasing.
T20 Remark 1 (projection under plan-anchored promises): without waiting the schedule from the mean times is p + D,
    so pi = W (D - beta)^+ under ETA promises and pi = 0 for D <= beta under day-definite promises.

Deviations are Gaussian, so g_r has the closed form sum_i w_i sigma_i psi((x - m_i)/sigma_i),
psi(z) = z Phi(z) + phi(z), m_i = s_i - mu_i.

Usage:
    python -I code/resequencing/theory_checks.py
"""
from __future__ import annotations

import itertools
import os
import sys

import numpy as np
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from avg_bound import T, iso_resid, quad_bound, tail_bound, tail_bound_lp  # noqa: E402

RNG = np.random.default_rng(2031)


def psi(z):
    return z * norm.cdf(z) + norm.pdf(z)


def g(x, m, sd, w):
    x = np.atleast_1d(x)[:, None]
    return np.sum(w * sd * psi((x - m) / sd), axis=1)


def G(y, m, sd, w):
    y = np.atleast_1d(y)[:, None]
    return np.sum(w * norm.cdf((y - m) / sd), axis=1)


def random_route(n=None, window=True):
    n = n or int(RNG.integers(3, 9))
    s = RNG.uniform(-30, 120, n) if window else np.full(n, 60.0)
    mu = RNG.normal(0, 5, n)
    sd = RNG.uniform(5, 40, n)
    return s - mu, sd, np.ones(n)


def eps_star(D, rho):
    """Direct O(N^2) evaluation of 1/2 max over pairs with D_a <= D_b of (rho_a - rho_b)^+."""
    best = 0.0
    for a in range(len(D)):
        for b in range(len(D)):
            if D[a] <= D[b]:
                best = max(best, rho[a] - rho[b])
    return best / 2


def topk_val(rho, score, K):
    return rho[np.argsort(-score, kind="stable")[:K]].sum()


def t1():
    worst = 0.0
    grid = np.linspace(-400, 400, 160001)
    for _ in range(50):
        m, sd, w = random_route()
        Gy = G(grid, m, sd, w)
        integ = np.concatenate([[0], np.cumsum((Gy[1:] + Gy[:-1]) / 2 * np.diff(grid))])
        xs = RNG.uniform(-100, 200, 20)
        worst = max(worst, float(np.max(np.abs(np.interp(xs, grid, integ) - g(xs, m, sd, w)))))
    print(f"T1 representation: max |integrated CDF - closed form| = {worst:.2e} -> {'PASS' if worst < 1e-3 else 'FAIL'}")
    return worst < 1e-3


def t2():
    ok_i = ok_w = True
    grid = np.linspace(-500, 600, 22001)
    for _ in range(400):
        N = int(RNG.integers(2, 12))
        K = int(RNG.integers(1, N))
        n = int(RNG.integers(3, 8))  # equal remaining counts -> equal total weight W = n
        routes = [random_route(n, window=bool(RNG.random() < 0.7)) for _ in range(N)]
        D = RNG.normal(0, 30, N)
        rho = np.array([g(D[r], *routes[r])[0] for r in range(N)])
        loss = topk_val(rho, rho, K) - topk_val(rho, D, K)
        e = eps_star(D, rho)
        ok_i &= loss <= 2 * min(K, N - K) * e + 1e-9
        Gs = np.array([G(grid, *rt) for rt in routes])
        Gbar = Gs.mean(axis=0)
        w1max = max(float(np.trapezoid(np.abs(Gs[r] - Gbar), grid)) for r in range(N))  # = W * W1(nu_r, nu_bar)
        ok_w &= loss <= 2 * min(K, N - K) * w1max + 1e-6
    print(f"T2 Proposition 2 and Corollary 1: {'PASS' if ok_i else 'FAIL'}; Wasserstein form: {'PASS' if ok_w else 'FAIL'}")
    return ok_i and ok_w


def t3():
    ok, ratios, constr = True, [], []
    xs = np.linspace(-120, 160, 141)
    for _ in range(60):
        N = int(RNG.integers(2, 4))
        n = int(RNG.integers(3, 7))
        routes = [random_route(n) for _ in range(N)]
        gx = np.array([g(xs, *rt) for rt in routes])
        Gamma = float(np.max(gx.max(axis=0) - gx.min(axis=0)))
        best = 0.0
        grid_idx = range(0, len(xs), 4)
        for combo in itertools.product(grid_idx, repeat=N):
            D = xs[list(combo)] + 1e-6 * np.arange(N)  # tiny offsets break ties deterministically
            rho = gx[np.arange(N), list(combo)]
            best = max(best, rho.max() - rho[np.argmax(D)])
        ok &= best <= Gamma + 1e-9
        ratios.append(best / Gamma if Gamma > 0 else 1.0)
        # constructive lower bound from the proof: r' at x, r just above x, all other routes far behind
        j = int(np.argmax(gx.max(axis=0) - gx.min(axis=0)))
        rp, r = int(np.argmax(gx[:, j])), int(np.argmin(gx[:, j]))
        D = np.full(N, -1e4)
        D[rp], D[r] = xs[j], xs[j] + 1e-6
        rho = np.array([g(D[k], *routes[k])[0] for k in range(N)])
        constr.append((rho.max() - rho[np.argmax(D)]) / Gamma if Gamma > 0 else 1.0)
    print(f"T3 Proposition 3, worst case K = 1: brute force never exceeds Gamma ({'PASS' if ok else 'FAIL'}); "
          f"coarse-grid brute force / Gamma: min {min(ratios):.3f}, median {np.median(ratios):.3f}; "
          f"proof construction / Gamma: min {min(constr):.6f}")
    return ok and min(constr) > 0.999


def t4():
    worst = -np.inf
    for _ in range(300):
        n = int(RNG.integers(1, 8))
        S = 20000
        D = RNG.normal(0, 30)
        s = RNG.uniform(-30, 120, n)
        xi = RNG.gamma(2.0, RNG.uniform(5, 20), size=(S, n)) - RNG.uniform(0, 30, n)  # skewed, any law
        mu_hat = xi.mean(axis=0) + RNG.normal(0, 10, n)  # any deterministic projection, not only the mean
        A, Ah = D + xi, D + mu_hat
        rho = np.mean(np.sum(np.maximum(A - s, 0), axis=1))
        pi = np.sum(np.maximum(Ah - s, 0))
        tau = np.mean(np.sum(np.abs(A - Ah) * (np.maximum(A, Ah) > s), axis=1))
        worst = max(worst, abs(rho - pi) - tau)
    print(f"T4 Proposition 5, route-level step, any projection: max(|rho - pi| - tau) = {worst:.2e} -> {'PASS' if worst <= 1e-9 else 'FAIL'}")
    return worst <= 1e-9


def t5():
    ok_chain = ok_dom = ok_lp = ok_iso = True
    for it in range(1000):
        N = int(RNG.integers(2, 60))
        K = int(RNG.integers(1, N))
        m = min(K, N - K)
        D = np.round(RNG.normal(0, 30, N), 0 if it % 3 == 0 else 3)  # every third instance has many delay ties
        rho = np.maximum(0.0, RNG.gamma(0.6, 20, N) + 0.3 * np.maximum(D, 0) * RNG.uniform(0, 2, N))
        loss = topk_val(rho, rho, K) - topk_val(rho, D, K)
        e_iso = iso_resid(D, rho)
        # a random nondecreasing gbar evaluated at D (piecewise constant on sorted distinct delays)
        u, grp = np.unique(D, return_inverse=True)
        gbar = np.cumsum(RNG.exponential(3, len(u))) + RNG.normal(0, 10)
        for e in (e_iso, rho - gbar[grp]):
            tb, qb = tail_bound(e, m), quad_bound(e, m)
            ok_chain &= (loss <= tb + 1e-9) and (tb <= qb + 1e-9)
            ok_dom &= tb <= 2 * m * np.max(np.abs(e)) + 1e-9
        ok_iso &= np.linalg.norm(e_iso) <= np.linalg.norm(rho - gbar[grp]) + 1e-9
        if it < 300:
            lp, _ = tail_bound_lp(D, rho, m)
            ok_lp &= (loss <= lp + 1e-6) and (lp <= tail_bound(e_iso, m) + 1e-6)
    ok = ok_chain and ok_dom and ok_lp and ok_iso
    print(f"T5 Proposition 4: chain {'PASS' if ok_chain else 'FAIL'}; below Prop 2 {'PASS' if ok_dom else 'FAIL'}; "
          f"LP between loss and isotonic tail {'PASS' if ok_lp else 'FAIL'}; isotonic minimizes ||e|| {'PASS' if ok_iso else 'FAIL'}")
    return ok


def t6():
    worst = 0.0
    for N, K in [(10, 3), (10, 5), (10, 7), (40, 2), (41, 30), (7, 6)]:
        m = min(K, N - K)
        c = 4.0
        D = np.concatenate([np.full(N - 2 * m - max(0, 2 * K - N), -1.0), np.zeros(m), np.ones(m),
                            np.full(max(0, 2 * K - N), 2.0)])
        rho = np.concatenate([np.zeros(N - 2 * m - max(0, 2 * K - N)), np.full(m, 2 * c), np.zeros(m),
                              np.full(max(0, 2 * K - N), 2 * c)])
        loss = topk_val(rho, rho, K) - topk_val(rho, D, K)
        e = iso_resid(D, rho)
        tb, qb = tail_bound(e, m), quad_bound(e, m)
        worst = max(worst, abs(loss - 2 * c * m), abs(tb - loss), abs(qb - loss))
    print(f"T6 Proposition 4 tightness: max |bound - loss| = {worst:.2e} -> {'PASS' if worst < 1e-9 else 'FAIL'}")
    return worst < 1e-9


def t7():
    ok_pt = ok_b = True
    grid = np.linspace(-500, 600, 22001)
    for _ in range(300):
        N = int(RNG.integers(2, 25))
        K = int(RNG.integers(1, N))
        m = min(K, N - K)
        n = int(RNG.integers(3, 8))
        routes = [random_route(n, window=bool(RNG.random() < 0.7)) for _ in range(N)]
        D = RNG.normal(0, 30, N)
        rho = np.array([g(D[r], *routes[r])[0] for r in range(N)])
        gbar = np.mean([g(D, *rt) for rt in routes], axis=0)  # fleet-average exposure function at each D_r
        Gs = np.array([G(grid, *rt) for rt in routes])
        Wd = np.array([np.trapezoid(np.abs(Gs[r] - Gs.mean(axis=0)), grid) for r in range(N)])  # W * d_r
        ok_pt &= bool(np.all(np.abs(rho - gbar) <= Wd + 1e-6))
        loss = topk_val(rho, rho, K) - topk_val(rho, D, K)
        ok_b &= (loss <= T(Wd, 2 * m) + 1e-6) and (T(Wd, 2 * m) <= np.sqrt(2 * m) * np.linalg.norm(Wd) + 1e-9)
    print(f"T7 Corollary 2: pointwise |g_r - gbar| <= W d_r {'PASS' if ok_pt else 'FAIL'}; "
          f"loss <= W T_2m(d) <= W sqrt(2m)||d|| {'PASS' if ok_b else 'FAIL'}")
    return ok_pt and ok_b


def lg_model(t_done, c, q, phi, V, S=0.0, Sf=None):
    """Linear-Gaussian route-shock model: returns Var D, beta, Var(theta|D), lambda_i, Var U_i."""
    Lam, Q = t_done.sum(), (t_done ** 2).sum()
    varD = V * (Lam ** 2 * phi + Q * (1 - phi)) + S
    beta = Lam * phi * V / varD
    vth = phi * V * (Q * (1 - phi) * V + S) / varD
    Sf = np.zeros_like(c) if Sf is None else Sf
    return varD, beta, vth, 1 + c * beta, c ** 2 * vth + q * (1 - phi) * V + Sf


def t8():
    """Monte Carlo of the leg-level model vs the closed form of Lemma 3 (risk as a function of delay)."""
    worst_mc, worst_int = 0.0, 0.0
    for _ in range(20):
        Ld, nf = int(RNG.integers(3, 12)), int(RNG.integers(2, 6))
        t_done = RNG.uniform(60, 400, Ld)
        t_fut = RNG.uniform(60, 400, (nf, 3))  # three future legs to each promise, cumulative below
        c = np.cumsum(t_fut.sum(axis=1))
        q = np.cumsum((t_fut ** 2).sum(axis=1))
        phi, V = RNG.uniform(0.05, 0.95), RNG.uniform(0.05, 0.3)
        s = RNG.uniform(-200, 900, nf)
        w = RNG.uniform(0.5, 2, nf)
        varD, beta, vth, lam, vU = lg_model(t_done, c, q, phi, V)
        # Monte Carlo: draw theta and leg noise, compute D and future deviations; estimate E[lateness | D near x]
        M = 400000
        th = RNG.normal(0, np.sqrt(phi * V), M)
        D = th * t_done.sum() + RNG.normal(0, 1, M) * np.sqrt((1 - phi) * V * (t_done ** 2).sum())
        fut = th[:, None] * c[None, :] + RNG.normal(0, 1, (M, nf)) * np.sqrt((1 - phi) * V * q)[None, :]
        late = np.maximum(D[:, None] + fut - s[None, :], 0) @ w
        # regression of lateness on D within a narrow window around x
        x = float(RNG.normal(0, np.sqrt(varD)))
        sel = np.abs(D - x) < 0.02 * np.sqrt(varD)
        mc = late[sel].mean()
        se = late[sel].std() / np.sqrt(sel.sum()) + 0.01  # 0.01 s floor: near-zero risk gives se = 0
        cf = float(np.sum(w * np.sqrt(vU) * psi((lam * x - s) / np.sqrt(vU))))
        worst_mc = max(worst_mc, abs(mc - cf) / se)  # in Monte Carlo standard errors
        # integral of G-tilde (from far in the lower tail up to x) equals the closed form
        sdY = np.sqrt(vU) / lam
        grid = np.linspace(min(x, (s / lam).min()) - 12 * sdY.max(), x, 200001)
        Gt = np.sum(w * lam * norm.cdf((grid[:, None] - s / lam) / sdY), axis=1)
        integ = np.trapezoid(Gt, grid)
        worst_int = max(worst_int, abs(integ - cf) / max(cf, 1.0))
    ok = worst_mc < 4.5 and worst_int < 1e-3
    print(f"T8 Lemma 3 (linear learning): Monte Carlo vs closed form max |diff| = {worst_mc:.2f} s.e.; "
          f"integral of G-tilde vs closed form {worst_int:.1e} -> {'PASS' if ok else 'FAIL'}")
    return ok


def t9():
    ok, tight = True, []
    grid = np.linspace(-800, 1500, 23001)
    for _ in range(300):
        N = int(RNG.integers(3, 30))
        K = int(RNG.integers(1, N))
        m = min(K, N - K)
        n = int(RNG.integers(2, 7))
        lam = [RNG.uniform(1.0, 2.5, n) for _ in range(N)]
        w = [1.0 / l for l in lam]  # equal delay sensitivity: sum w * lambda = n for every route
        s = [RNG.uniform(-30, 300, n) for _ in range(N)]
        sdU = [RNG.uniform(5, 60, n) for _ in range(N)]
        D = RNG.normal(60, 80, N)
        rho = np.array([np.sum(w[r] * sdU[r] * psi((lam[r] * D[r] - s[r]) / sdU[r])) for r in range(N)])
        Wt = float(n)
        Gt = np.array([np.sum(w[r] * lam[r] * norm.cdf((grid[:, None] - s[r] / lam[r]) * lam[r] / sdU[r]), axis=1)
                       for r in range(N)])
        mu = np.array([np.sum(w[r] * lam[r] * (s[r] / lam[r])) / Wt for r in range(N)])
        dt = np.array([np.trapezoid(np.abs(Gt[r] - Gt.mean(axis=0)), grid) / Wt for r in range(N)])
        SD = np.argsort(-D, kind="stable")[:K]
        EK = float(np.mean(np.maximum(D[SD] - mu[SD], 0)))
        if EK <= 0:
            continue
        loss_rel = (topk_val(rho, rho, K) - topk_val(rho, D, K)) / topk_val(rho, rho, K)
        b1 = T(dt, 2 * m) / (K * EK)
        b2 = 2 * m / K * dt.max() / EK
        ok &= (loss_rel <= b1 + 1e-6) and (b1 <= b2 + 1e-12)
        tight.append(loss_rel / b1 if b1 > 0 else 0)
    print(f"T9 Proposition 6: bound holds {'PASS' if ok else 'FAIL'} (median loss/bound {np.median(tight):.3f})")
    return ok


def t10():
    ok = True
    for _ in range(200):
        Ld = int(RNG.integers(2, 60))
        t_done = RNG.uniform(30, 600, Ld)
        S = float(RNG.uniform(0, 2e4))
        V = RNG.uniform(0.05, 0.3)
        phis = np.linspace(0, 1, 51)
        out = [lg_model(t_done, np.array([1000.0]), np.array([2e5]), p, V, S) for p in phis]
        b = np.array([o[1] for o in out])
        vd = np.array([o[0] for o in out])
        ok &= bool(np.all(np.diff(b) > 0) and np.all(np.diff(vd) > 0))
        r0 = lg_model(t_done, np.array([1.0]), np.array([1.0]), 0.0, V)[0]
        r1 = lg_model(t_done, np.array([1.0]), np.array([1.0]), 1.0, V)[0]
        ok &= abs(r1 / r0 - t_done.sum() ** 2 / (t_done ** 2).sum()) < 1e-9
    # exact lognormal variance and Var(theta | D) against Monte Carlo
    worst_v, worst_th = 0.0, 0.0
    for _ in range(6):
        Ld = int(RNG.integers(3, 25))
        t_done = RNG.uniform(60, 500, Ld)
        phi, V = RNG.uniform(0.05, 0.95), RNG.uniform(0.05, 0.2)
        M = 400000
        th = RNG.normal(0, np.sqrt(phi * V), M)
        eps = RNG.normal(0, np.sqrt((1 - phi) * V), (M, Ld))
        D = (t_done[None, :] * (np.exp(th[:, None] + eps - V / 2) - 1)).sum(axis=1)
        Lam, Q = t_done.sum(), (t_done ** 2).sum()
        exact = (np.exp(phi * V) - 1) * (Lam ** 2 - Q) + (np.exp(V) - 1) * Q
        worst_v = max(worst_v, abs(D.var() / exact - 1))
        # linear-Gaussian Var(theta | D): residual variance of theta after linear regression on D
        Dl = Lam * th + (t_done[None, :] * eps).sum(axis=1)
        _, beta, vth, _, _ = lg_model(t_done, np.array([1.0]), np.array([1.0]), phi, V)
        worst_th = max(worst_th, abs(np.var(th - beta * Dl) / vth - 1))
    ok2 = worst_v < 0.02 and worst_th < 0.02
    print(f"T10 Proposition 7: monotone beta and Var D, ratio n_eff {'PASS' if ok else 'FAIL'}; exact lognormal Var D "
          f"(max rel. err {worst_v:.3f}) and Var(theta|D) (max rel. err {worst_th:.3f}) {'PASS' if ok2 else 'FAIL'}")
    return ok and ok2


def t11():
    worst = -np.inf
    for _ in range(20000):
        N = int(RNG.integers(2, 30))
        K = int(RNG.integers(1, N))
        m = min(K, N - K)
        rho = RNG.exponential(1, N) * (RNG.random(N) < 0.6)
        sig = rho + RNG.normal(0, RNG.uniform(0.01, 2), N)  # any score, e.g. a projection pi
        loss = np.sort(rho)[::-1][:K].sum() - topk_val(rho, sig, K)
        e = rho - sig
        tail = T(np.maximum(e, 0), m) + T(np.maximum(-e, 0), m)
        tau = np.abs(e) + RNG.exponential(0.1, N)  # any tau_r >= |rho_r - pi_r|
        worst = max(worst, loss - tail, tail - T(np.abs(e), 2 * m), T(np.abs(e), 2 * m) - 2 * m * np.abs(e).max(),
                    loss - T(tau, 2 * m))
    ok = worst <= 1e-9
    print(f"T11 Lemma 2 tail chain for arbitrary scores and Proposition 5 tail form: max violation {worst:.2e} -> "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def t12():
    ok = True
    for N in range(2, 25):
        for K in range(1, N):
            m = min(K, N - K)
            if K <= N - K:
                D = np.r_[np.ones(K), np.zeros(K), -np.ones(N - 2 * K)]
                rho = np.r_[np.zeros(K), np.ones(K), np.zeros(N - 2 * K)]
            else:
                D = np.r_[np.ones(N - K), np.zeros(N - K), 2 * np.ones(2 * K - N)]
                rho = np.r_[np.zeros(N - K), np.ones(N - K), np.ones(2 * K - N)]
            loss = np.sort(rho)[::-1][:K].sum() - topk_val(rho, D, K)
            ok &= abs(loss - 2 * m * eps_star(D, rho)) < 1e-12
    print(f"T12 Proposition 2 attained for all 2 <= N <= 24, 1 <= K < N: {'PASS' if ok else 'FAIL'}")
    return ok


def t13():
    ok = True
    lo, hi = -20.0, 40.0
    for _ in range(300):
        N = int(RNG.integers(3, 10))
        K = int(RNG.integers(1, N))
        routes = [random_route() for _ in range(N)]
        grid = np.linspace(lo + 1.0, hi - 1.0, 601)  # interior points of X = [lo, hi]
        vals = np.array([g(grid, *rt) for rt in routes])  # N x grid
        gap = vals.max(0) - vals.min(0)
        j = int(np.argmax(gap))
        x, rp, r = grid[j], int(np.argmax(vals[:, j])), int(np.argmin(vals[:, j]))
        if gap[j] <= 0:
            continue
        # Proposition 3 construction: T + {r} at x, everyone else (incl. r') at x - delta
        delta = 1e-7
        others = [k for k in range(N) if k not in (r, rp)]
        Tset = others[:K - 1]
        D = np.full(N, x - delta)
        D[Tset + [r]] = x
        rho = np.array([g(D[k], *routes[k])[0] for k in range(N)])
        SD = set(np.argsort(-D, kind="stable")[:K].tolist())
        loss = np.sort(rho)[::-1][:K].sum() - rho[list(SD)].sum()
        ok &= SD == set(Tset + [r]) and loss >= gap[j] - 1e-5 and loss <= min(K, N - K) * gap.max() + 1e-9
        # Proposition 1 construction: T at x0 + 2 eps, r at x0 + eps, others at x0; strictly suboptimal
        eps_ = 1e-3
        D1 = np.full(N, x)
        D1[Tset] = x + 2 * eps_
        D1[r] = x + eps_
        rho1 = np.array([g(D1[k], *routes[k])[0] for k in range(N)])
        SD1 = np.argsort(-D1, kind="stable")[:K]
        ok &= set(SD1.tolist()) == set(Tset + [r]) and np.sort(rho1)[::-1][:K].sum() > rho1[SD1].sum()
    print(f"T13 Propositions 1 and 3 on a bounded delay interval, any K: {'PASS' if ok else 'FAIL'}")
    return ok


def t14():
    x = RNG.exponential(3, 10**6)
    ok = all(bool(np.all(np.maximum(x - t, 0) <= x**2 / (4 * t) + 1e-12)) for t in (0.1, 1.0, 5.0, 20.0))

    def tbar(X, alpha):
        q = np.quantile(X, 1 - alpha)
        return alpha * q + np.maximum(X - q, 0).mean()

    ok &= all(tbar(x, a) <= np.sqrt(a * np.mean(x**2)) + 1e-9 for a in (0.05, 0.1, 0.25, 0.5))
    worst = -np.inf
    for kappa in (0.05, 0.2, 0.5, 0.8):
        N = 200000
        K = int(kappa * N)
        D = RNG.normal(0, 1, N)
        rho = np.exp(D) * RNG.lognormal(-0.5 * 0.6**2, 0.6, N) * (1 + 0.5 * (D < 0))
        gfix = np.exp(D)  # fixed nondecreasing g, not fitted
        e = rho - gfix
        alpha = min(kappa, 1 - kappa)
        loss_per = (np.sort(rho)[::-1][:K].sum() - topk_val(rho, D, K)) / N
        bound = tbar(np.maximum(e, 0), alpha) + tbar(np.maximum(-e, 0), alpha)
        ok &= loss_per <= bound * 1.01 and bound <= np.sqrt(2 * alpha * np.mean(e**2)) + 1e-9
        worst = max(worst, loss_per / bound)
    print(f"T14 Corollary 3: elementary inequality, Tbar bound, large-fleet limit (max loss/bound {worst:.3f}) -> "
          f"{'PASS' if ok else 'FAIL'}")
    return bool(ok)


def t15():
    ok, worst = True, 0.0
    for s, txt in ((0.25, 0.199), (0.50, 0.395), (1.00, 0.766)):
        X = np.exp(s * RNG.standard_normal(4_000_000) - s ** 2 / 2)
        mc = float(np.mean(np.abs(X - 1)))
        cf = 2 * (2 * norm.cdf(s / 2) - 1)
        worst = max(worst, abs(mc - cf))
        ok &= abs(mc - cf) < 3e-3 and abs(cf - txt) < 5e-4 and abs(X.mean() - 1) < 3e-3
    print(f"T15 Lemma 5: E|X - 1| closed form vs Monte Carlo max |diff| = {worst:.1e} -> {'PASS' if ok else 'FAIL'}")
    return bool(ok)


def t16():
    from gate3_calibrate import re_decomp
    Jr = RNG.integers(2, 40, 60)  # unbalanced legs per route
    J, G = Jr.sum(), len(Jr)
    J0 = (J - (Jr ** 2).sum() / J) / (G - 1)
    pr, pl = 0.3, 1.1  # psi_route^2, psi_leg^2
    msw, msb = [], []
    grp = np.repeat(np.arange(G), Jr)
    for _ in range(4000):
        z = 0.7 + np.repeat(RNG.normal(0, np.sqrt(pr), G), Jr) + RNG.normal(0, np.sqrt(pl), J)
        zr = np.bincount(grp, z) / Jr
        msw.append(np.sum((z - zr[grp]) ** 2) / (J - G))
        msb.append(np.sum(Jr * (zr - z.mean()) ** 2) / (G - 1))
    ok = abs(np.mean(msw) / pl - 1) < 0.01 and abs(np.mean(msb) / (pl + J0 * pr) - 1) < 0.03
    out = re_decomp(z, grp)
    ok &= abs(out["kappa2"] - msw[-1]) < 1e-9 and abs(out["tau2"] - max(0.0, (msb[-1] - msw[-1]) / J0)) < 1e-9
    print(f"T16 Lemma 6: E MSW / psi_leg^2 = {np.mean(msw) / pl:.3f}, E MSB / (psi_leg^2 + J0 psi_route^2) = "
          f"{np.mean(msb) / (pl + J0 * pr):.3f}; matches re_decomp -> {'PASS' if ok else 'FAIL'}")
    return bool(ok)


def t17():
    """Proposition 5 tightness (appendix proof): group 1 has Ahat - d = c, A = d; group 2 has
    Ahat - d = c - eps, A - d = 2c - eps; the rest have Ahat = A (late by 3c if K > N - K, early otherwise)."""
    ok, c, eps = True, 1.0, 1e-3
    for N in range(2, 25):
        for K in range(1, N):
            m = min(K, N - K)
            rest = N - 2 * m
            fill = 3 * c if K > N - K else -1.0
            ahat = np.r_[np.full(m, c), np.full(m, c - eps), np.full(rest, fill)]  # projected arrival - deadline
            a = np.r_[np.zeros(m), np.full(m, 2 * c - eps), np.full(rest, fill)]  # actual arrival - deadline
            rho, pi = np.maximum(a, 0), np.maximum(ahat, 0)
            tau = np.abs(a - ahat) * (np.maximum(a, ahat) > 0)
            srt = np.sort(pi)[::-1]
            ok &= K == N or srt[K - 1] > srt[K]  # no ties at the selection boundary
            ok &= bool(np.all(np.abs(rho - pi) <= tau + 1e-12))
            loss = np.sort(rho)[::-1][:K].sum() - topk_val(rho, pi, K)
            t2m = np.sort(tau)[::-1][:2 * m].sum()
            ok &= abs(t2m - 2 * m * tau.max()) < 1e-12 and loss <= t2m + 1e-12
            ok &= abs(loss - m * (2 * c - eps)) < 1e-9 and abs(loss / t2m - (1 - eps / (2 * c))) < 1e-12
    print(f"T17 Proposition 5 tight for all 2 <= N <= 24, 1 <= K < N (loss / bound = 1 - eps/2c): "
          f"{'PASS' if ok else 'FAIL'}")
    return bool(ok)


def schedule(t0, T, S, o, wait):
    """Recursion (9): T[..., k] travel into remaining stop k, S[..., k] service there; waits for o[k] if wait[k]."""
    A = np.empty(np.broadcast(T, S).shape)
    cur = t0 + T[..., 0]
    A[..., 0] = np.where(wait[0], np.maximum(cur, o[0]), cur)
    for k in range(1, T.shape[-1]):
        cur = A[..., k - 1] + S[..., k - 1] + T[..., k]
        A[..., k] = np.where(wait[k], np.maximum(cur, o[k]), cur)
    return A


def t18():
    ok, worst_jensen, ratios = True, -np.inf, []
    for _ in range(40):
        N, n, M = int(RNG.integers(6, 30)), int(RNG.integers(3, 10)), 400
        rho, pi, tau = np.empty(N), np.empty(N), np.empty(N)
        for r in range(N):
            tbar, sbar = RNG.uniform(5, 30, n), RNG.uniform(1, 8, n)
            T = tbar * np.exp(RNG.normal(0, 0.2, (M, 1)) + RNG.normal(0, 0.3, (M, n)))
            S = sbar * np.exp(RNG.normal(0, 0.4, (M, n)))
            o, wait = RNG.uniform(0, 150, n), RNG.random(n) < 0.5
            d = o + RNG.uniform(0, 60, n)
            A = schedule(0.0, T, S, o, wait)  # empirical law: M equally likely scenarios
            Ahat = schedule(0.0, T.mean(0), S.mean(0), o, wait)  # schedule from the (empirical) mean times
            worst_jensen = max(worst_jensen, float(np.max(Ahat - A.mean(0))))
            rho[r] = np.maximum(A - d, 0).mean(0).sum()
            pi[r] = np.maximum(Ahat - d, 0).sum()
            tau[r] = (np.abs(A - Ahat) * (np.maximum(A, Ahat) > d)).mean(0).sum()
        ok &= bool(np.all(pi <= rho + 1e-9))
        for K in range(1, N):
            m = min(K, N - K)
            loss = np.sort(rho)[::-1][:K].sum() - topk_val(rho, pi, K)
            bound = np.sort(tau)[::-1][:m].sum()
            ok &= loss <= bound + 1e-9
            if bound > 0:
                ratios.append(loss / bound)
    ok &= worst_jensen <= 1e-9
    c, eps = 1.0, 1e-3
    for N in range(2, 25):
        for K in range(1, N):
            m = min(K, N - K)
            rest = N - 2 * m
            fill = 3 * c if K > N - K else -1.0
            ahat = np.r_[np.full(m, c), np.full(m, c - eps), np.full(rest, fill)]
            a = np.r_[np.full(m, c), np.full(m, 2 * c - eps), np.full(rest, fill)]
            rho, pi = np.maximum(a, 0), np.maximum(ahat, 0)
            tau = np.abs(a - ahat) * (np.maximum(a, ahat) > 0)
            srt = np.sort(pi)[::-1]
            ok &= bool(np.all(ahat <= a)) and srt[K - 1] > srt[K]
            loss = np.sort(rho)[::-1][:K].sum() - topk_val(rho, pi, K)
            tm = np.sort(tau)[::-1][:m].sum()
            ok &= abs(tm - m * tau.max()) < 1e-12 and loss <= tm + 1e-12 and abs(loss / tm - (1 - eps / c)) < 1e-12
    print(f"T18 Corollary 4: max(Ahat - E A) = {worst_jensen:.1e}; loss <= T_m(tau) on random fleets (max loss/bound "
          f"{max(ratios):.3f}); construction tight -> {'PASS' if ok else 'FAIL'}")
    return bool(ok)


def loss_fn(u):
    """Standard normal loss function L(u) = phi(u) - u (1 - Phi(u)) = E (Z - u)^+."""
    return norm.pdf(u) - u * norm.sf(u)


def t19():
    """Lemma 4: dispersion term and projection error of one promise when A - Ahat ~ N(0, sigma^2)."""
    ok, worst = True, 0.0
    Z = RNG.standard_normal(4_000_000)
    for z in (-2.5, -1.0, -0.3, 0.0, 0.4, 1.0, 2.0, 3.0):
        sigma = float(RNG.uniform(5, 30))
        ahat, d = 100.0, 100.0 + z * sigma  # projected slack of z standard deviations
        A = ahat + sigma * Z
        tau_mc = float(np.mean(np.abs(A - ahat) * (np.maximum(A, ahat) > d)))
        err_mc = float(np.mean(np.maximum(A - d, 0)) - max(ahat - d, 0.0))
        tau_cf = sigma * norm.pdf(z) if z >= 0 else sigma * np.sqrt(2 / np.pi)
        err_cf = sigma * loss_fn(abs(z))
        worst = max(worst, abs(tau_mc - tau_cf) / sigma, abs(err_mc - err_cf) / sigma)
        ok &= abs(tau_mc - tau_cf) < 2e-3 * sigma and abs(err_mc - err_cf) < 2e-3 * sigma
        ok &= err_cf <= tau_cf + 1e-12 and (z >= 0 or err_cf <= sigma * norm.pdf(0) + 1e-12)
    u = np.linspace(0, 6, 601)
    ok &= bool(np.all(np.diff(loss_fn(u)) < 0)) and bool(np.all((loss_fn(u) >= 0) & (loss_fn(u) <= norm.pdf(u))))
    print(f"T19 Lemma 4: Gaussian dispersion term and projection error, closed form vs Monte Carlo max |diff|/sigma = "
          f"{worst:.1e}; L decreasing, 0 <= L <= phi -> {'PASS' if ok else 'FAIL'}")
    return bool(ok)


def t20():
    """Remark 1: without waiting, the schedule from the mean times at check-in time P0 + D is p + D, so under ETA
    promises (d = p + beta) pi = W (D - beta)^+, and under day-definite promises (d = P + beta) pi = 0 if D <= beta."""
    ok, beta = True, 60.0
    for _ in range(200):
        n = int(RNG.integers(3, 40))
        tbar, sbar = RNG.uniform(2, 25, n), RNG.uniform(1, 6, n)
        P0 = float(RNG.uniform(0, 300))  # planned end of service at the stop just served
        no_wait = np.zeros(n, bool)
        p = schedule(P0, tbar, sbar, np.zeros(n), no_wait)  # planned arrivals at the remaining stops
        P_end = p[-1] + sbar[-1]  # planned route end
        w = RNG.integers(0, 2, n).astype(float)
        for D in RNG.uniform(-90, 180, 5):
            ahat = schedule(P0 + D, tbar, sbar, np.zeros(n), no_wait)
            ok &= bool(np.allclose(ahat, p + D))
            pi_eta = float(np.sum(w * np.maximum(ahat - (p + beta), 0)))
            ok &= abs(pi_eta - w.sum() * max(D - beta, 0.0)) < 1e-8
            pi_day = float(np.sum(w * np.maximum(ahat - (P_end + beta), 0)))
            ok &= D > beta or pi_day == 0.0
    print(f"T20 Remark 1: projection under ETA promises = W (D - beta)^+, under day-definite promises 0 for D <= beta "
          f"-> {'PASS' if ok else 'FAIL'}")
    return bool(ok)


if __name__ == "__main__":
    res = [t1(), t2(), t3(), t4(), t5(), t6(), t7(), t8(), t9(), t10(), t11(), t12(), t13(), t14(), t15(), t16(), t17(), t18(),
           t19(), t20()]
    print("ALL PASS" if all(res) else "SOME CHECK FAILED")
