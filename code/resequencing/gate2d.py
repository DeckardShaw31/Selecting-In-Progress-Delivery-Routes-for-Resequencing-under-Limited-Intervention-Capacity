"""Gate 2d: observed-lateness check of behind-schedule vs. promise-aware triage
(Konovalenko, Hvattum & Iversen 2026 data). Implements analysis/gate2d/preregistration.md.

Usage:
    python -I code/resequencing/gate2d.py <stops.pkl> <out_dir>
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

FRACS = (0.35, 0.50)
KS = (0.05, 0.10)
SHRINK = 10.0
N_BOOT, BOOT_SEED = 2000, 2028


def fit_pace(d, routes):
    """OLS of consecutive actual inter-arrival times on DistanceA over all legs of the given routes
    (actual order by IndexA), keeping legs with 0 < dt < 99th percentile."""
    x = d[d["Route ID"].isin(routes)].sort_values(["Route ID", "IndexA"])
    dt = x.groupby("Route ID")["Arrived Time"].diff()
    dist = x["DistanceA"]
    m = dt.notna() & (dt > 0)
    q = dt[m].quantile(0.99)
    m &= dt < q
    X = np.column_stack([np.ones(m.sum()), dist[m].to_numpy()])
    y = dt[m].to_numpy()
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    r2 = 1 - np.sum((y - X @ coef) ** 2) / np.sum((y - y.mean()) ** 2)
    return float(coef[0]), float(coef[1]), float(r2), int(m.sum())


def project(rem, t0, a, b):
    """Planned-order projection of the remaining stops: each leg a + b*DistanceP; wait until Earliest."""
    t = t0
    late = 0.0
    for _, s in rem.iterrows():
        t += a + b * s["DistanceP"]
        if not np.isnan(s["Earliest Time"]) and t < s["Earliest Time"]:
            t = s["Earliest Time"]
        late += max(0.0, t - s["Latest Time"])
    return late


def route_scores(g_all, f, a, b):
    """Scores and outcome for one route at check-in fraction f."""
    g_all = g_all.sort_values("IndexA")
    nd = g_all[g_all["Depot"] == 0]
    n = len(nd)
    m = int(round(f * n))
    if m < 1 or n - m < 1:
        return None
    served = nd.iloc[:m]
    remaining = nd.iloc[m:]
    lat_served = (served["Arrived Time"] - served["Latest Time"]).clip(lower=0)
    lat_rem = (remaining["Arrived Time"] - remaining["Latest Time"]).clip(lower=0)
    t0 = float(served["Arrived Time"].iloc[-1])
    rem_planned = remaining.sort_values("IndexP")
    # route-specific intercept from completed legs (all rows in actual order up to the m-th served stop)
    last_pos = served["IndexA"].iloc[-1]
    done = g_all[g_all["IndexA"] <= last_pos]
    dt = done["Arrived Time"].diff().iloc[1:]
    dist = done["DistanceA"].iloc[1:]
    ok = dt > 0
    a_hat = float(np.mean(dt[ok] - b * dist[ok])) if ok.any() else a
    a_r = (m * a_hat + SHRINK * a) / (m + SHRINK)
    return dict(
        n_stops=n, m=m, f=f,
        B1_late_so_far=float(lat_served.sum()),
        B2_currently_late=float(lat_served.iloc[-1]),
        F_projection=project(rem_planned, t0, a, b),
        F_learn=project(rem_planned, t0, a_r, b),
        Y_future_late_min=float(lat_rem.sum()),
        Y_future_late_count=int((lat_rem > 0).sum()),
        median_rem_window=float((remaining["Latest Time"] - remaining["Earliest Time"]).median()),
    )


def value(df, score, k):
    K = max(1, int(round(k * len(df))))
    order = np.argsort(-df[score].to_numpy(), kind="stable")
    return df["Y_future_late_min"].to_numpy()[order[:K]].sum()


def main():
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    d = pd.read_pickle(src)
    nd = d[d["Depot"] == 0]
    per = nd.groupby("Route ID").size()
    incomplete = nd.groupby("Route ID").apply(lambda g: g["Arrived Time"].isna().any() or g["Latest Time"].isna().any())
    elig = set(per[per >= 10].index) - set(incomplete[incomplete].index)
    route_week = d.groupby("Route ID")["Week ID"].first()
    weeks = sorted(d["Week ID"].unique())
    calib_weeks = set(weeks[: len(weeks) // 2])
    calib = [r for r in elig if route_week[r] in calib_weeks]
    evalr = sorted(r for r in elig if route_week[r] not in calib_weeks)
    a, b, r2, nlegs = fit_pace(d, calib)
    info = dict(eligible=len(elig), calibration_routes=len(calib), evaluation_routes=len(evalr),
                calibration_weeks=len(calib_weeks), pace_a_min=a, pace_b_min_per_km=b, pace_r2=r2, pace_legs=nlegs)
    print(json.dumps(info, indent=1))
    rows = []
    groups = d[d["Route ID"].isin(evalr)].groupby("Route ID")
    for rid, g in groups:
        for f in FRACS:
            s = route_scores(g, f, a, b)
            if s:
                s["route"] = rid
                rows.append(s)
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(out, "gate2d_routes.csv"), index=False)
    # efficiencies
    tab = []
    for f in FRACS:
        x = res[res.f == f].sort_values("route").reset_index(drop=True)
        for k in KS:
            vo = value(x, "Y_future_late_min", k)
            rec = dict(f=f, k=k, routes=len(x), V_oracle=vo)
            for sc in ["B1_late_so_far", "B2_currently_late", "F_projection", "F_learn"]:
                rec[sc] = value(x, sc, k) / vo if vo > 0 else np.nan
            rec["random"] = x.Y_future_late_min.mean() * max(1, int(round(k * len(x)))) / vo
            rec["ratio_B1_over_F"] = rec["B1_late_so_far"] / rec["F_projection"]
            tab.append(rec)
    tab = pd.DataFrame(tab)
    tab.to_csv(os.path.join(out, "gate2d_efficiency.csv"), index=False)
    print(tab.round(3).to_string(index=False))
    # D1 with bootstrap
    x = res[res.f == 0.35].sort_values("route").reset_index(drop=True)
    rng = np.random.default_rng(BOOT_SEED)
    diffs = []
    for _ in range(N_BOOT):
        bs = x.iloc[rng.integers(0, len(x), len(x))].reset_index(drop=True)
        vo = value(bs, "Y_future_late_min", 0.05)
        diffs.append((value(bs, "F_projection", 0.05) - value(bs, "B1_late_so_far", 0.05)) / vo if vo > 0 else np.nan)
    lo, hi = np.nanpercentile(diffs, [2.5, 97.5])
    t = tab[(tab.f == 0.35) & (tab.k == 0.05)].iloc[0]
    ratio = t.ratio_B1_over_F
    d1 = (ratio <= 0.85) and (lo > 0)
    outcome = "SUPPORT" if d1 else ("CONSISTENT-WEAK" if ratio <= 1 else "CONTRADICT")
    print(f"\nD1: f=0.35, k=5%: V(B1)/V(F) = {ratio:.3f} (<= 0.85); bootstrap 95% CI of [V(F)-V(B1)]/V(oracle) = [{lo:.3f}, {hi:.3f}] (> 0)")
    print(f"Outcome: {outcome}")
    # exploratory dose analogue: terciles of median remaining window width
    x["tercile"] = pd.qcut(x.median_rem_window.rank(method="first"), 3, labels=["narrow", "middle", "wide"])
    ter = []
    for tname, g in x.groupby("tercile", observed=True):
        g = g.reset_index(drop=True)
        vo = value(g, "Y_future_late_min", 0.05)
        ter.append(dict(tercile=tname, routes=len(g), median_window=g.median_rem_window.median(),
                        B1=value(g, "B1_late_so_far", 0.05) / vo, F=value(g, "F_projection", 0.05) / vo,
                        ratio_B1_over_F=value(g, "B1_late_so_far", 0.05) / max(value(g, "F_projection", 0.05), 1e-9)))
    ter = pd.DataFrame(ter)
    ter.to_csv(os.path.join(out, "gate2d_terciles.csv"), index=False)
    print("\nExploratory (window-width terciles, f=0.35, k=5%):")
    print(ter.round(3).to_string(index=False))
    json.dump(dict(info, D1_ratio=ratio, D1_ci=[lo, hi], outcome=outcome), open(os.path.join(out, "gate2d_summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
