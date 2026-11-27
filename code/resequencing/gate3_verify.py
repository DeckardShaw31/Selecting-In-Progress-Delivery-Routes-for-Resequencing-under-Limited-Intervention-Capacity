"""Gate 3 verification checks 1-3 (check 4, the Proposition 2 bounds, runs inside analyze_gate3.py).

1. With sigma_U = 0, the base split and seed 2027, gate3.run_config reproduces the Gate 2c rows
   (gate1.py output) on 20 routes, both noise levels, both check-ins, all 7 regimes.
2. Lateness computed from start_times equals eval_seq on random instances; and the route-level
   Proposition 2(ii) inequality |L2 - rho| <= tau holds on the check-1 rows.
3. The closed form for eps* matches a brute-force L-infinity isotonic fit on random data with ties.

Usage:
    python -I code/resequencing/gate3_verify.py <train_routes.pkl> <gate2c_main.csv> [--workers W]
"""
from __future__ import annotations

import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gate3  # noqa: E402
from analyze_gate3 import eps_star  # noqa: E402
from gate1 import NOISE  # noqa: E402
from sim import eval_seq  # noqa: E402

COMPARE = ["rho_min", "lam_min", "delta_min", "L2_proj", "L3_prior", "delay_plan_min", "b_late_stops",
           "b_late_window_stops", "stops_moved", "passes", "rho_in_min"]


def check1(routes_pkl, g2c_csv, workers):
    g2c = pd.read_csv(g2c_csv)
    ridx = sorted(g2c.ridx.unique())[:: max(1, len(g2c.ridx.unique()) // 20)][:20]
    tasks = []
    for r in ridx:
        for nz in ("base", "high"):
            sR, sE, sS = NOISE[nz]
            tasks.append(dict(ridx=int(r), seed=2027, arm="A", noise=nz, rng_level=gate3.NOISE_IDX[nz], sR=sR, sE=sE,
                              sS=sS, phi=np.nan, sU=0.0, fracs=[0.35, 0.50], regimes=gate3.REG_W, misspec=False))
    with Pool(workers, initializer=gate3._init, initargs=(routes_pkl,)) as pool:
        rows = sum(pool.map(gate3.run_config, tasks), [])
    new = pd.DataFrame(rows)
    old = g2c[g2c.ridx.isin(ridx)]
    m = new.merge(old, on=["ridx", "frac", "regime", "noise"], suffixes=("", "_g2c"))
    assert (m.route == m.route_g2c).all()
    print(f"check 1: {len(new)} new rows, {len(old)} Gate 2c rows, {len(m)} matched")
    worst = 0.0
    for c in COMPARE:
        dmax = float(np.max(np.abs(m[c].astype(float) - m[f"{c}_g2c"].astype(float))))
        worst = max(worst, dmax)
        print(f"   max |diff| {c:22s} {dmax:.3e}")
    ok1 = len(m) == len(old) == len(new) and worst < 1e-6
    ineq = np.abs(new.L2_proj - new.rho_min) <= new.tau_ii_min + 1e-6
    print(f"check 1 {'PASS' if ok1 else 'FAIL'}; route-level Prop 2(ii) inequality holds in {ineq.mean():.1%} of rows "
          f"({'PASS' if ineq.all() else 'FAIL'})")
    return ok1 and bool(ineq.all())


def check2(n_inst=300, seed=11):
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(n_inst):
        N = int(rng.integers(3, 12))
        S = int(rng.integers(1, 20))
        T = rng.uniform(100, 900, size=(S, N, N))
        Sv = rng.uniform(0, 300, size=(S, N))
        Sv[:, 0] = 0
        opn = np.where(rng.random(N) < 0.4, rng.uniform(0, 4000, N), -1e18)
        cls = np.where(rng.random(N) < 0.7, rng.uniform(500, 6000, N), 1e18)
        w = rng.uniform(0, 2, N)
        w[0] = 0
        wait = bool(rng.random() < 0.5)
        seq = (rng.permutation(N - 1) + 1).astype(np.int64)
        t0 = float(rng.uniform(0, 1000))
        A = gate3.start_times(seq, t0, T, Sv, opn, wait)
        late = np.mean(np.sum(w[seq][None, :] * np.maximum(A - cls[seq][None, :], 0), axis=1))
        ref = eval_seq(seq, t0, T, Sv, opn, cls, w, wait, 0.0)
        worst = max(worst, abs(late - ref) / max(1.0, abs(ref)))
    ok = worst < 1e-9
    print(f"check 2: {n_inst} random instances, max relative |diff| = {worst:.2e} -> {'PASS' if ok else 'FAIL'}")
    return ok


def brute_eps(D, rho, iters=80):
    """Smallest eps with a nondecreasing g, |rho - g(D)| <= eps (bisection with a greedy feasibility test)."""
    gb = pd.DataFrame({"D": D, "r": rho}).groupby("D")["r"].agg(["max", "min"]).sort_index()
    hi_v, lo_v = gb["max"].to_numpy(), gb["min"].to_numpy()

    def feasible(e):
        g = -np.inf
        for mx, mn in zip(hi_v, lo_v):
            g = max(g, mx - e)  # lowest admissible value for this group
            if g > mn + e:
                return False
        return True

    lo, hi = 0.0, float(rho.max() - rho.min()) + 1.0
    for _ in range(iters):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if feasible(mid) else (mid, hi)
    return hi


def check3(n_inst=300, seed=12):
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(n_inst):
        n = int(rng.integers(2, 60))
        D = rng.integers(0, max(2, n // 2), size=n).astype(float)  # ties on purpose
        rho = rng.exponential(10, size=n) * (rng.random(n) < 0.7)
        worst = max(worst, abs(eps_star(D, rho) - brute_eps(D, rho)))
    ok = worst < 1e-6
    print(f"check 3: {n_inst} random instances, max |closed form - brute force| = {worst:.2e} -> {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    w = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 14
    ok2 = check2()
    ok3 = check3()
    ok1 = check1(sys.argv[1], sys.argv[2], w)
    print("ALL PASS" if (ok1 and ok2 and ok3) else "SOME CHECK FAILED")
