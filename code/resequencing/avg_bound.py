"""Average-case loss bounds for behind-schedule triage (Proposition 5 in analysis/theory/proof_audit.md, Section 3b; Proposition 4 in the manuscript).

For risks rho, delays D, capacity K, m = min(K, N - K), any nondecreasing gbar and residuals e = rho - gbar(D):

    L = V(S_rho) - V(S_D) <= T_m(e+) + T_m(e-) <= sqrt(2 m) ||e||_2,

where T_m(x) is the sum of the m largest entries of x >= 0. With gbar the least-squares isotonic regression of rho
on D, ||e||_2^2 = N s_rho^2 H_iso. The tail bound minimized over gbar is a linear program (tail_bound_lp).

The evaluation on simulated cells follows analysis/theory/average_bound_check.md, which was written before this
script was first run on any cell.

Usage:
    python -I code/resequencing/avg_bound.py <gate2c_main.csv> <gate3_A.csv> <gate3_B.csv> <out_dir>
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import linprog
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression

SHARE = 0.05


def topk(target, score, K):
    """Sum of target over the top-K entries of score (stable: ties keep route order), as in analyze_gate3.py."""
    return float(target[np.argsort(-score, kind="stable")[:K]].sum())


def T(x, m):
    """Sum of the m largest entries of a nonnegative vector."""
    if m <= 0:
        return 0.0
    return float(np.sort(x)[::-1][:m].sum())


def iso_resid(D, rho):
    fit = IsotonicRegression(increasing=True, out_of_bounds="clip").fit(D, rho).predict(D)
    return rho - fit


def tail_bound(e, m):
    return T(np.maximum(e, 0), m) + T(np.maximum(-e, 0), m)


def quad_bound(e, m):
    return float(np.sqrt(2 * m) * np.linalg.norm(e))


def eps_star(D, rho):
    gb = pd.DataFrame({"D": D, "r": rho}).groupby("D")["r"].agg(["max", "min"]).sort_index()
    pm = np.maximum.accumulate(gb["max"].to_numpy())
    return 0.5 * float(np.max(pm - gb["min"].to_numpy()))


def tail_bound_lp(D, rho, m):
    """min over nondecreasing gbar of T_m(e+) + T_m(e-), using T_m(x+) = min_{t >= 0} m t + sum (x_i - t)+.

    Variables: f (one per distinct delay, nondecreasing), u, v (N each, >= 0), t1, t2 (>= 0).
    Returns (value, residuals at the optimum).
    """
    N = len(rho)
    uniq, grp = np.unique(D, return_inverse=True)
    U = len(uniq)
    nv = U + 2 * N + 2
    iu, iv, it1, it2 = U, U + N, U + 2 * N, U + 2 * N + 1
    c = np.zeros(nv)
    c[iu:iu + N] = 1.0
    c[iv:iv + N] = 1.0
    c[it1] = c[it2] = m
    r = np.arange(N)
    # rho_i - f_g(i) - t1 <= u_i   ->  -f_g(i) - t1 - u_i <= -rho_i
    A1 = sparse.csr_matrix((np.concatenate([-np.ones(N), -np.ones(N), -np.ones(N)]),
                            (np.concatenate([r, r, r]), np.concatenate([grp, np.full(N, it1), iu + r]))), shape=(N, nv))
    # f_g(i) - rho_i - t2 <= v_i   ->   f_g(i) - t2 - v_i <= rho_i
    A2 = sparse.csr_matrix((np.concatenate([np.ones(N), -np.ones(N), -np.ones(N)]),
                            (np.concatenate([r, r, r]), np.concatenate([grp, np.full(N, it2), iv + r]))), shape=(N, nv))
    # f_k - f_{k+1} <= 0
    k = np.arange(U - 1)
    A3 = sparse.csr_matrix((np.concatenate([np.ones(U - 1), -np.ones(U - 1)]),
                            (np.concatenate([k, k]), np.concatenate([k, k + 1]))), shape=(max(U - 1, 0), nv))
    A = sparse.vstack([A1, A2, A3]).tocsr()
    b = np.concatenate([-rho, rho, np.zeros(max(U - 1, 0))])
    bounds = [(None, None)] * U + [(0, None)] * (2 * N + 2)
    res = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method="highs")
    if not res.success:
        raise RuntimeError(res.message)
    f = res.x[:U]
    return float(res.fun), rho - f[grp]


def cell_stats(D, rho, share=SHARE, lp=True):
    N = len(rho)
    K = max(1, int(round(share * N)))
    m = min(K, N - K)
    v_opt = topk(rho, rho, K)
    loss = v_opt - topk(rho, D, K)
    e = iso_resid(D, rho)
    H = float(np.sum(e ** 2) / np.sum((rho - rho.mean()) ** 2))
    out = dict(N=N, K=K, V_opt=v_opt, loss_rel=loss / v_opt, H_iso=H,
               B_inf=2 * m * eps_star(D, rho) / v_opt,
               B_tail_iso=tail_bound(e, m) / v_opt,
               B_quad_iso=quad_bound(e, m) / v_opt,
               sd_over_topmean=float(rho.std() / (v_opt / K)))
    if lp:
        out["B_tail_opt"] = tail_bound_lp(D, rho, m)[0] / v_opt
    tol = 1e-9
    out["holds"] = bool(out["loss_rel"] <= min(out["B_tail_iso"], out.get("B_tail_opt", np.inf), out["B_inf"]) + tol
                        and out["B_tail_iso"] <= out["B_quad_iso"] + tol)
    return out


def cells(df, keys):
    rows = []
    for key, g in df.groupby(keys):
        g = g.sort_values("ridx")
        rec = dict(zip(keys, key if isinstance(key, tuple) else (key,)))
        rec.update(cell_stats(g.delay_plan_min.to_numpy(float), g.rho_min.to_numpy(float)))
        rows.append(rec)
    return pd.DataFrame(rows)


def summarize(t, name):
    out = dict(set=name, cells=len(t), all_hold=bool(t.holds.all()))
    for b in ("B_quad_iso", "B_tail_iso", "B_tail_opt", "B_inf", "H_iso"):
        out[f"spearman_{b}"] = float(spearmanr(t[b], t.loss_rel).statistic)
        if b.startswith("B_"):
            out[f"nonvacuous_{b}"] = float((t[b] < 1).mean())
            out[f"median_ratio_{b}"] = float(np.median(t[b] / t.loss_rel.clip(lower=1e-12)))
    return out


def main(gate2c, gate3a, gate3b, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    cols = ["ridx", "delay_plan_min", "rho_min"]
    c2 = cells(pd.read_csv(gate2c, usecols=cols + ["regime", "frac", "noise"]), ["regime", "frac", "noise"])
    a = pd.read_csv(gate3a, usecols=cols + ["sU", "regime", "frac", "noise"])
    aw = cells(a[a.sU == 0], ["regime", "frac", "noise"])
    au = cells(a[a.sU > 0], ["sU", "regime", "noise"])
    b = cells(pd.read_csv(gate3b, usecols=cols + ["noise", "phi", "regime"]), ["noise", "phi", "regime"])
    sets = {"gate2c (pre-specified)": c2, "gate3_armB (pre-specified)": b,
            "gate3_W (exploratory)": aw, "gate3_sigmaU (exploratory)": au}
    for name, t in sets.items():
        t.to_csv(os.path.join(out_dir, f"avg_bound_{name.split()[0]}.csv"), index=False)
    summ = [summarize(t, n) for n, t in sets.items()]
    json.dump(summ, open(os.path.join(out_dir, "avg_bound_summary.json"), "w"), indent=1)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        for name, t in sets.items():
            print(f"\n{name}")
            print(t.drop(columns=["V_opt", "holds"]).round(3).to_string(index=False))
        print()
        print(pd.DataFrame(summ).set_index("set").T.to_string())


if __name__ == "__main__":
    main(*sys.argv[1:5])
