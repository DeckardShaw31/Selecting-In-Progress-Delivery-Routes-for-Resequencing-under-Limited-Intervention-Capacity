"""Gate 2d follow-ups (do not change the pre-registered decision):
  (a) pre-registered SECONDARY outcome: count of remaining stops arriving after Latest;
  (b) EXPLORATORY robustness (post hoc, logged): exclude routes with any non-depot stop late by > 24 h
      (timestamp errors: max arrival ~12.8 million minutes; 653 visits late by > 24 h).

Usage:
    python -I code/resequencing/gate2d_robustness.py <stops.pkl> <gate2d_out_dir>
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

SCORES = ["B1_late_so_far", "B2_currently_late", "F_projection", "F_learn"]


def value(df, score, k, ycol):
    K = max(1, int(round(k * len(df))))
    order = np.argsort(-df[score].to_numpy(), kind="stable")
    return df[ycol].to_numpy()[order[:K]].sum()


def table(res, ycol, label):
    rows = []
    for f in (0.35, 0.50):
        x = res[res.f == f].sort_values("route").reset_index(drop=True)
        for k in (0.05, 0.10):
            vo = value(x, ycol, k, ycol)
            rec = dict(set=label, outcome=ycol, f=f, k=k, routes=len(x))
            for sc in SCORES:
                rec[sc] = value(x, sc, k, ycol) / vo
            rec["random"] = x[ycol].mean() * max(1, int(round(k * len(x)))) / vo
            rec["ratio_B1_over_F"] = rec["B1_late_so_far"] / rec["F_projection"]
            rows.append(rec)
    return pd.DataFrame(rows)


def boot(res, ycol, f=0.35, k=0.05, n=2000, seed=2028):
    x = res[res.f == f].sort_values("route").reset_index(drop=True)
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        bs = x.iloc[rng.integers(0, len(x), len(x))].reset_index(drop=True)
        vo = value(bs, ycol, k, ycol)
        out.append((value(bs, "F_projection", k, ycol) - value(bs, "B1_late_so_far", k, ycol)) / vo if vo > 0 else np.nan)
    return [float(v) for v in np.nanpercentile(out, [2.5, 97.5])]


def main():
    src, out = sys.argv[1], sys.argv[2]
    res = pd.read_csv(os.path.join(out, "gate2d_routes.csv"))
    d = pd.read_pickle(src)
    nd = d[d["Depot"] == 0]
    bad = set(nd.loc[(nd["Arrived Time"] - nd["Latest Time"]) > 1440, "Route ID"].unique())
    clean = res[~res.route.isin(bad)]
    tabs = [table(res, "Y_future_late_count", "all routes (pre-registered secondary)"),
            table(clean, "Y_future_late_min", "excl. >24h-late routes (exploratory)"),
            table(clean, "Y_future_late_count", "excl. >24h-late routes (exploratory)")]
    t = pd.concat(tabs, ignore_index=True)
    t.to_csv(os.path.join(out, "gate2d_followups.csv"), index=False)
    with pd.option_context("display.width", 250):
        print(t.round(3).to_string(index=False))
    ci = dict(secondary_count_all=boot(res, "Y_future_late_count"),
              exploratory_minutes_clean=boot(clean, "Y_future_late_min"),
              exploratory_count_clean=boot(clean, "Y_future_late_count"),
              excluded_routes=len(bad & set(res.route)), evaluation_routes=int(res.route.nunique()))
    print(json.dumps(ci, indent=1))
    json.dump(ci, open(os.path.join(out, "gate2d_followups_ci.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
