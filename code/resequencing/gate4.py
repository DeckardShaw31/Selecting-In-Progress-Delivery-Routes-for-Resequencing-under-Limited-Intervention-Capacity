"""Gate 4: second observed-outcome check on LaDe pickup tasks. Implements analysis/gate4/preregistration.md.

Usage:
    python -I code/resequencing/gate4.py <lade_pickup.pkl> <out_dir>
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
from numba import njit

FRACS = (0.35, 0.50)
KS = (0.05, 0.10)
SHRINK = 10.0
N_BOOT, BOOT_SEED = 2000, 2030
SCORES = ["B1_late_so_far", "B1c_late_count_so_far", "B2_currently_late", "F_projection", "F_nowait", "F_direct", "F_learn"]


def hav(lat1, lng1, lat2, lng2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi, dl = p2 - p1, np.radians(lng2 - lng1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * 6371.0088 * np.arcsin(np.sqrt(np.minimum(1.0, a)))


@njit(cache=True)
def _hav(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lo2 - lo1) / 2) ** 2
    return 2 * 6371.0088 * np.arcsin(np.sqrt(min(1.0, a)))


@njit(cache=True)
def project(t0, la0, lo0, la, lo, ws, we, a, b, wait):
    """Projected late minutes of the tasks in the given order from (la0, lo0) at t0."""
    t, pla, plo, late = t0, la0, lo0, 0.0
    for i in range(la.shape[0]):
        t += a + b * _hav(pla, plo, la[i], lo[i])
        if wait and t < ws[i]:
            t = ws[i]
        if t > we[i]:
            late += t - we[i]
        pla, plo = la[i], lo[i]
    return late


def route_bounds(df):
    k = (df.city.astype(str) + "|" + df.courier_id.astype(str) + "|" + df.ds.astype(str)).to_numpy()
    starts = np.flatnonzero(np.r_[True, k[1:] != k[:-1]])
    ends = np.r_[starts[1:], len(df)]
    return starts, ends


def fit_pace(df):
    """Per city: OLS of consecutive inter-pickup times on the distance between consecutive tasks,
    calibration courier-days, legs with 0 < dt < city 99th percentile."""
    out = {}
    for city, g in df.groupby("city", sort=True):
        s, e = route_bounds(g)
        pk, la, lo = g.pickup_time.to_numpy(), g.lat.to_numpy(), g.lng.to_numpy()
        same = np.ones(len(g), bool)
        same[s] = False  # first task of each route has no incoming leg
        dt = np.r_[np.nan, np.diff(pk)]
        dist = np.r_[np.nan, hav(la[:-1], lo[:-1], la[1:], lo[1:])]
        m = same & (dt > 0)
        q = np.quantile(dt[m], 0.99)
        m &= dt < q
        X = np.column_stack([np.ones(m.sum()), dist[m]])
        coef, *_ = np.linalg.lstsq(X, dt[m], rcond=None)
        r2 = 1 - np.sum((dt[m] - X @ coef) ** 2) / np.sum((dt[m] - dt[m].mean()) ** 2)
        out[city] = dict(a=float(coef[0]), b=float(coef[1]), r2=float(r2), legs=int(m.sum()), dt_q99=float(q))
    return out


def score_routes(df, pace):
    s, e = route_bounds(df)
    city = df.city.to_numpy()
    cour, ds = df.courier_id.to_numpy(), df.ds.to_numpy()
    acc, ws, we, pk = (df[c].to_numpy() for c in ("accept_time", "time_window_start", "time_window_end", "pickup_time"))
    la, lo = df.lat.to_numpy(), df.lng.to_numpy()
    rows = []
    for r0, r1 in zip(s, e):
        n = r1 - r0
        a, b = pace[city[r0]]["a"], pace[city[r0]]["b"]
        late_all = np.maximum(pk[r0:r1] - we[r0:r1], 0.0)
        for f in FRACS:
            m = int(round(f * n))
            if m < 1 or n - m < 1:
                continue
            t0 = pk[r0 + m - 1]
            la0, lo0 = la[r0 + m - 1], lo[r0 + m - 1]
            rem = np.arange(r0 + m, r1)
            rem = rem[acc[rem] <= t0]
            if len(rem) == 0:
                continue
            d0 = hav(la0, lo0, la[rem], lo[rem])
            order = rem[np.lexsort((d0, we[rem]))]
            served = late_all[:m]
            # route-specific intercept from completed legs among the m served tasks
            if m >= 2:
                dt = np.diff(pk[r0:r0 + m])
                dist = hav(la[r0:r0 + m - 1], lo[r0:r0 + m - 1], la[r0 + 1:r0 + m], lo[r0 + 1:r0 + m])
                ok = dt > 0
                a_hat = float(np.mean(dt[ok] - b * dist[ok])) if ok.any() else a
            else:
                a_hat = a
            a_r = (m * a_hat + SHRINK * a) / (m + SHRINK)
            yl = np.maximum(pk[rem] - we[rem], 0.0)
            rows.append(dict(
                city=city[r0], courier_id=cour[r0], ds=ds[r0], n_tasks=n, m=m, f=f, known_remaining=len(rem),
                B1_late_so_far=float(served.sum()), B1c_late_count_so_far=int((served > 0).sum()),
                B2_currently_late=float(served[-1]),
                F_projection=project(t0, la0, lo0, la[order], lo[order], ws[order], we[order], a, b, True),
                F_nowait=project(t0, la0, lo0, la[order], lo[order], ws[order], we[order], a, b, False),
                F_direct=float(np.sum(np.maximum(t0 + a + b * d0 - we[rem], 0.0))),
                F_learn=project(t0, la0, lo0, la[order], lo[order], ws[order], we[order], a_r, b, True),
                Y_count=int((yl > 0).sum()), Y_minutes=float(yl.sum()),
                median_rem_window=float(np.median(we[rem] - ws[rem])),
                unknown_share=float(np.mean(acc[r0:r1] > t0)),
            ))
    return pd.DataFrame(rows)


def value(df, score, k, y):
    K = max(1, int(round(k * len(df))))
    return df[y].to_numpy()[np.argsort(-df[score].to_numpy(), kind="stable")[:K]].sum()


def eff_rows(x, y, label):
    out = []
    for k in KS:
        vo = value(x, y, k, y)
        rec = dict(subset=label, outcome=y, f=x.f.iloc[0], k=k, routes=len(x), V_oracle=vo)
        for sc in SCORES:
            rec[sc] = value(x, sc, k, y) / vo if vo > 0 else np.nan
        rec["random"] = x[y].mean() * max(1, int(round(k * len(x)))) / vo if vo > 0 else np.nan
        rec["ratio_B1_over_F"] = rec["B1_late_so_far"] / rec["F_projection"] if rec["F_projection"] > 0 else np.nan
        out.append(rec)
    return out


def boot(x, y, seed=BOOT_SEED):
    rng = np.random.default_rng(seed)
    diffs = []
    for _ in range(N_BOOT):
        bs = x.iloc[rng.integers(0, len(x), len(x))].reset_index(drop=True)
        vo = value(bs, y, 0.05, y)
        diffs.append((value(bs, "F_projection", 0.05, y) - value(bs, "B1_late_so_far", 0.05, y)) / vo if vo > 0 else np.nan)
    return [float(v) for v in np.nanpercentile(diffs, [2.5, 97.5])]


def label(ratio, lo):
    return "SUPPORT" if (ratio <= 0.85 and lo > 0) else ("CONSISTENT-WEAK" if ratio <= 1 else "CONTRADICT")


def main():
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    d = pd.read_pickle(src)
    pace = fit_pace(d[d.split == "calibration"])
    print("pace per city:", json.dumps(pace, indent=1))
    ev = d[d.split == "evaluation"]
    res = score_routes(ev, pace)
    res = res.sort_values(["city", "courier_id", "ds", "f"]).reset_index(drop=True)
    res.to_csv(os.path.join(out, "gate4_routes.csv"), index=False)
    tab = []
    for f in FRACS:
        x = res[res.f == f].reset_index(drop=True)
        for y in ("Y_count", "Y_minutes"):
            tab += eff_rows(x, y, "all")
            for c, g in x.groupby("city"):
                tab += eff_rows(g.reset_index(drop=True), y, f"city={c}")
    tab = pd.DataFrame(tab)
    tab.to_csv(os.path.join(out, "gate4_efficiency.csv"), index=False)
    show = ["subset", "outcome", "f", "k", "routes"] + SCORES + ["random", "ratio_B1_over_F"]
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(tab[tab.subset == "all"][show].round(3).to_string(index=False))
    x = res[res.f == 0.35].reset_index(drop=True)
    summ = dict(pace=pace, evaluation_routes_f035=len(x), evaluation_routes_f050=int((res.f == 0.5).sum()))
    for tag, y in (("D1", "Y_count"), ("D2", "Y_minutes")):
        t = tab[(tab.subset == "all") & (tab.outcome == y) & (tab.f == 0.35) & (tab.k == 0.05)].iloc[0]
        lo, hi = boot(x, y)
        summ[tag] = dict(outcome=y, ratio=float(t.ratio_B1_over_F), ci=[lo, hi], label=label(t.ratio_B1_over_F, lo))
        print(f"{tag} ({y}): f=0.35, k=5%: V(B1)/V(F) = {t.ratio_B1_over_F:.3f} (<= 0.85); 95% CI of [V(F)-V(B1)]/V(oracle) = [{lo:.3f}, {hi:.3f}] -> {summ[tag]['label']}")
    summ["outcome"] = summ["D1"]["label"]
    print(f"Outcome (D1): {summ['outcome']}")
    json.dump(summ, open(os.path.join(out, "gate4_summary.json"), "w"), indent=1)
    # ---- reported without thresholds (after the decision is recorded)
    wd = []
    for y in ("Y_count", "Y_minutes"):
        num = {sc: 0.0 for sc in SCORES}
        den = 0.0
        for _, g in x.groupby(["city", "ds"]):
            K = max(1, int(round(0.05 * len(g))))
            yy = g[y].to_numpy()
            den += yy[np.argsort(-yy, kind="stable")[:K]].sum()
            for sc in SCORES:
                num[sc] += yy[np.argsort(-g[sc].to_numpy(), kind="stable")[:K]].sum()
        wd.append(dict(outcome=y, **{sc: num[sc] / den for sc in SCORES}))
    pd.DataFrame(wd).to_csv(os.path.join(out, "gate4_within_cityday.csv"), index=False)
    ter = []
    for col, name in (("median_rem_window", "window_width"), ("unknown_share", "unknown_share")):
        x["ter"] = pd.qcut(x[col].rank(method="first"), 3, labels=["low", "middle", "high"])
        for tn, g in x.groupby("ter", observed=True):
            g = g.reset_index(drop=True)
            vb, vf = value(g, "B1_late_so_far", 0.05, "Y_count"), value(g, "F_projection", 0.05, "Y_count")
            ter.append(dict(split_by=name, tercile=tn, routes=len(g), median_value=float(g[col].median()),
                            ratio_B1_over_F=vb / vf if vf > 0 else np.nan))
    pd.DataFrame(ter).to_csv(os.path.join(out, "gate4_terciles.csv"), index=False)
    ev_tasks = ev
    desc = dict(share_tasks_late=float((ev_tasks.pickup_time > ev_tasks.time_window_end).mean()),
                share_tasks_early=float((ev_tasks.pickup_time < ev_tasks.time_window_start).mean()),
                window_width_min_median=float((ev_tasks.time_window_end - ev_tasks.time_window_start).median()),
                window_width_min_p5_p95=[float(v) for v in (ev_tasks.time_window_end - ev_tasks.time_window_start).quantile([0.05, 0.95])],
                tasks_per_courier_day_median=float(x.n_tasks.median()),
                known_remaining_median_f035=float(x.known_remaining.median()),
                unknown_share_median_f035=float(x.unknown_share.median()),
                routes_with_future_late_f035=float((x.Y_count > 0).mean()))
    json.dump(desc, open(os.path.join(out, "gate4_descriptives.json"), "w"), indent=1)
    print("\nwithin city-day:\n", pd.DataFrame(wd).round(3).to_string(index=False))
    print("\nterciles:\n", pd.DataFrame(ter).round(3).to_string(index=False))
    print("\ndescriptives:", json.dumps(desc, indent=1))


if __name__ == "__main__":
    main()
