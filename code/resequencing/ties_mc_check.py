"""Two implementation diagnostics reported in Section 4.3 (added in the fourth review pass, after the gates).

1. Tie-breaking. Every selection keeps the K highest scores with a stable sort, so tied scores are broken by the
   fixed order of the routes in the result files. For every cell and score of Gates 2c, 3 (arms A and B), 2d and 4,
   this script computes the exact expected efficiency when ties at the selection boundary are broken uniformly at
   random instead, and reports the largest change: over all cells, and over the efficiencies the paper reports.
2. Monte Carlo edge of L4. The L4 score is the 256-scenario held-out estimate of rho, which also values delta. Scoring
   L4 instead on the independent 48-scenario optimization estimate (rho_in) bounds the edge this sharing gives L4.

    python -I code/resequencing/ties_mc_check.py           (run from the project root; writes analysis/theory/ties_mc/)
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.join(ROOT, "analysis", "theory", "ties_mc")


def topk(y, s, K):
    return float(y[np.argsort(-s, kind="stable")[:K]].sum())


def random_ties(y, s, K):
    """Expected sum of y over the K highest s when ties at the K-th value are broken uniformly at random."""
    order = np.argsort(-s, kind="stable")
    thr = s[order[K - 1]]
    above, tie = s > thr, s == thr
    return float(y[above].sum() + (K - above.sum()) * y[tie].mean())


def tie_rows():
    rows = []

    def add(source, cell, y, scores, denom, shares, reported):
        for k in shares:
            K = max(1, int(round(k * len(y))))
            v = topk(y, denom, K)
            if v <= 0:
                continue
            for name, s in scores.items():
                s = np.asarray(s, float)
                if np.isnan(s).any():
                    continue
                a, b = topk(y, s, K) / v, random_ties(y, s, K) / v
                rows.append(dict(source=source, cell=cell, score=name, k=k, stable=a, random=b, diff=abs(a - b),
                                 reported=reported(cell, name, k)))

    # Gate 2c: Table 4 and the capacity sentence report L1-L3 at 2, 5 and 10%.
    d = pd.read_csv(os.path.join(ROOT, "analysis/gate2c/results/gate2c_main.csv"))
    for key, g in d.groupby(["regime", "noise", "frac"]):
        g = g.reset_index(drop=True)
        rep = (lambda cell, name, k: k == 0.05 or name == "L1"
               or (name == "L2" and cell[0] in ("MIXB", "MIXB20", "MIXB50", "TW") and cell[1] == "base"))
        add("gate2c", key, g.delta_min.to_numpy(float),
            {"L1": g.delay_plan_min, "L2": g.L2_proj, "L3": g.L3_prior}, g.rho_min.to_numpy(float),
            (0.02, 0.05, 0.10), rep)
    # Gate 3: the paper reports arm A at 5% (all ladder levels) and arm B for L1 and L3 at 5%.
    cols = ("delay_plan_min", "L2_proj", "L3_prior", "L4f", "L3_m05", "L4_m05", "L3_m2", "L4_m2")
    for arm in ("A", "B"):
        d = pd.read_csv(os.path.join(ROOT, f"analysis/gate3/results/gate3_{arm}.csv"))
        keys = [c for c in ("sU", "phi", "regime", "noise", "frac") if c in d.columns]
        for key, g in d.groupby(keys):
            g = g.reset_index(drop=True)
            rep = (lambda cell, name, k, arm=arm: k == 0.05
                   and (arm == "A" or name in ("delay_plan_min", "L3_prior")))
            add("gate3_" + arm, key, g.delta_min.to_numpy(float), {c: g[c] for c in cols if c in g.columns},
                g.rho_min.to_numpy(float), (0.02, 0.05, 0.10), rep)
    # Observed data: Table 5 reports B1, B2 and F at 5% and 10%, check-in 35%; the text also cites 50%.
    obs = ("B1_late_so_far", "B2_currently_late", "F_projection")
    g4 = pd.read_csv(os.path.join(ROOT, "analysis/gate4/gate4_routes.csv"))
    for f, g in g4.groupby("f"):
        g = g.reset_index(drop=True)
        for yc in ("Y_count", "Y_minutes"):
            y = g[yc].to_numpy(float)
            add("gate4", (f, yc), y, {c: g[c] for c in obs}, y, (0.05, 0.10), lambda cell, name, k: True)
    g2 = pd.read_csv(os.path.join(ROOT, "analysis/gate2d/gate2d_routes.csv"))
    for f, g in g2.groupby("f"):
        g = g.sort_values("route").reset_index(drop=True)  # the order gate2d.py uses
        for yc in ("Y_future_late_min", "Y_future_late_count"):
            y = g[yc].to_numpy(float)
            add("gate2d", (f, yc), y, {c: g[c] for c in obs}, y, (0.05, 0.10), lambda cell, name, k: True)
    return pd.DataFrame(rows)


def mc_rows():
    rows = []
    d = pd.read_csv(os.path.join(ROOT, "analysis/gate2c/results/gate2c_main.csv"))
    a = pd.read_csv(os.path.join(ROOT, "analysis/gate3/results/gate3_A.csv"))
    for source, frame, keys in (("gate2c", d, ["regime", "noise", "frac"]),
                                ("gate3_A_sU0", a[a.sU == 0], ["regime", "noise", "frac"])):
        for key, g in frame.groupby(keys):
            y, K = g.delta_min.to_numpy(float), max(1, round(0.05 * len(g)))
            rows.append(dict(source=source, cell=key,
                             ratio=topk(y, g.rho_in_min.to_numpy(float), K) / topk(y, g.rho_min.to_numpy(float), K)))
    return pd.DataFrame(rows)


def main():
    os.makedirs(OUT, exist_ok=True)
    t = tie_rows()
    m = mc_rows()
    t.to_csv(os.path.join(OUT, "ties.csv"), index=False)
    m.to_csv(os.path.join(OUT, "mc_edge.csv"), index=False)
    worst_all = t.sort_values("diff", ascending=False).iloc[0]
    worst_rep = t[t.reported].sort_values("diff", ascending=False).iloc[0]
    summ = dict(
        ties_cells_checked=int(len(t)),
        ties_max_change_all=float(worst_all["diff"]), ties_worst_all=str(worst_all[["source", "cell", "score", "k"]].tolist()),
        ties_max_change_reported=float(worst_rep["diff"]),
        ties_worst_reported=str(worst_rep[["source", "cell", "score", "k"]].tolist()),
        mc_cells=int(len(m)), mc_min_ratio=float(m.ratio.min()), mc_median_ratio=float(m.ratio.median()),
    )
    json.dump(summ, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))
    ok = summ["ties_max_change_reported"] <= 0.03 and summ["mc_min_ratio"] >= 0.98
    print("text claims hold (reported change <= 0.03; L4 on independent draws >= 0.98):", "PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
