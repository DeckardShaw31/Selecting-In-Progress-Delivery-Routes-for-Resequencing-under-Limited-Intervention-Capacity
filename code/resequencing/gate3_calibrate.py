"""Gate 3 calibration facts, measured before the Gate 3 pre-registration was frozen.

Descriptive only: no triage score or lateness outcome is computed here.

1. Station-day sizes in the Amazon splits (routes per station and date), to decide whether the
   Gate 3 sample can consist of whole station-days.
2. Dispersion of observed inter-arrival times in the Konovalenko et al. calibration weeks:
   one-way random-effects decomposition of log(actual inter-arrival / fleet-pace prediction)
   into a route effect (tau^2) and a leg residual (kappa^2).
3. The same decomposition on simulated executions of Amazon training routes at the assumed noise
   levels, with and without a persistent leg-level error in the carrier's mean times (sigma_U).

Usage:
    python -I code/resequencing/gate3_calibrate.py <project_root> <out_json>
"""
from __future__ import annotations

import json
import os
import pickle
import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate1 import NOISE  # noqa: E402
from gate2d import fit_pace  # noqa: E402


def station_days(route_data_path):
    rd = json.load(open(route_data_path, encoding="utf-8"))
    all_c, elig_c = Counter(), Counter()
    for v in rd.values():
        key = (v["station_code"], v["date_YYYY_MM_DD"])
        all_c[key] += 1
        if len(v["stops"]) - 1 >= 10:
            elig_c[key] += 1
    sizes = np.array(list(elig_c.values()))
    tot = sizes.sum()
    return dict(
        routes=int(sum(all_c.values())), eligible_routes=int(tot), station_days=len(elig_c),
        stations=len({k[0] for k in elig_c}), dates=len({k[1] for k in elig_c}),
        size_median=float(np.median(sizes)), size_p25=float(np.percentile(sizes, 25)),
        size_p75=float(np.percentile(sizes, 75)), size_max=int(sizes.max()),
        share_routes_in_days_ge10=float(sizes[sizes >= 10].sum() / tot),
        share_routes_in_days_ge20=float(sizes[sizes >= 20].sum() / tot),
        days_ge20=int((sizes >= 20).sum()),
    )


def re_decomp(resid, groups):
    """Method-of-moments one-way random-effects decomposition (unbalanced)."""
    df = pd.DataFrame({"r": resid, "g": groups})
    grp = df.groupby("g")["r"]
    n_i = grp.size().to_numpy().astype(float)
    m_i = grp.mean().to_numpy()
    N, G = n_i.sum(), len(n_i)
    grand = df.r.mean()
    ssw = float(((df.r - grp.transform("mean")) ** 2).sum())
    ssb = float((n_i * (m_i - grand) ** 2).sum())
    msw = ssw / (N - G)
    msb = ssb / (G - 1)
    n0 = (N - (n_i ** 2).sum() / N) / (G - 1)
    tau2 = max(0.0, (msb - msw) / n0)
    return dict(legs=int(N), routes=int(G), mean=float(grand), tau2=tau2, kappa2=msw,
                tau=float(np.sqrt(tau2)), kappa=float(np.sqrt(msw)),
                route_share=float(tau2 / (tau2 + msw)), total_sd=float(np.sqrt(tau2 + msw)))


def konovalenko(stops_path):
    d = pd.read_pickle(stops_path)
    nd = d[d["Depot"] == 0]
    per = nd.groupby("Route ID").size()
    incomplete = nd.groupby("Route ID").apply(lambda g: g["Arrived Time"].isna().any() or g["Latest Time"].isna().any())
    elig = set(per[per >= 10].index) - set(incomplete[incomplete].index)
    route_week = d.groupby("Route ID")["Week ID"].first()
    weeks = sorted(d["Week ID"].unique())
    calib_weeks = set(weeks[: len(weeks) // 2])
    calib = sorted(r for r in elig if route_week[r] in calib_weeks)
    a, b, r2, nlegs = fit_pace(d, calib)
    x = d[d["Route ID"].isin(calib)].sort_values(["Route ID", "IndexA"]).copy()
    x["dt"] = x.groupby("Route ID")["Arrived Time"].diff()
    nonmono = set(x.loc[x.dt < 0, "Route ID"])
    ok = x.dt.notna() & (x.dt > 0)
    q99 = x.loc[ok, "dt"].quantile(0.99)
    ok &= (x.dt < q99) & (x["Depot"] == 0)
    y = x[ok]
    resid = np.log(y.dt.to_numpy()) - np.log(a + b * y["DistanceA"].to_numpy())
    out = dict(pace_a=a, pace_b=b, pace_r2=r2, pace_legs=nlegs, calibration_routes=len(calib),
               dt_q99=float(q99), all=re_decomp(resid, y["Route ID"].to_numpy()))
    keep = ~y["Route ID"].isin(nonmono).to_numpy()
    out["monotone_routes_only"] = re_decomp(resid[keep], y["Route ID"].to_numpy()[keep])
    out["routes_nonmonotone"] = len(nonmono)
    return out


PHI_GRID = (0.10, 0.25, 0.50, 0.75, 0.90)


def simulated(train_pkl, n_routes=300, seed=3031):
    routes = pickle.load(open(train_pkl, "rb"))
    elig = [i for i, r in enumerate(routes) if r["n"] - 1 >= 10]
    pick = np.random.default_rng(seed).choice(elig, size=n_routes, replace=False)
    settings = [(f"{nz}_sU{sU:.2f}", NOISE[nz], sU) for nz in ("low", "base", "high") for sU in (0.0, 0.25, 0.50)]
    # Gate 3 shock-split grid: route-wide share phi of travel log-variance at the base/high total
    for lvl in ("base", "high"):
        sR0, sE0, sS0 = NOISE[lvl]
        V = sR0 ** 2 + sE0 ** 2
        for phi in PHI_GRID:
            settings.append((f"split_{lvl}_phi{phi:.2f}", (np.sqrt(phi * V), np.sqrt((1 - phi) * V), sS0), 0.0))
    res = {}
    for label, (sR, sE, sS), sU in settings:
        resid, grp = [], []
        for j, i in enumerate(pick):
            r = routes[i]
            rng = np.random.default_rng([seed, int(i), int(sU * 100), len(label)])
            n = r["n"]
            tb = np.array([float(r["tbar"][k - 1, k]) for k in range(1, n)])
            sb = r["sbar"][: n - 1].astype(float)  # service at the previous position (0 at the station)
            u = rng.normal(0.0, sU, size=n - 1) if sU > 0 else np.zeros(n - 1)
            v = rng.normal(0.0, sU, size=n - 1) if sU > 0 else np.zeros(n - 1)
            t_true = tb * np.exp(u - sU ** 2 / 2)
            s_true = sb * np.exp(v - sU ** 2 / 2)
            theta = rng.normal(0.0, sR)
            legs = t_true * np.exp(theta + rng.normal(0.0, sE, size=n - 1) - (sR ** 2 + sE ** 2) / 2)
            svc = s_true * np.exp(rng.normal(0.0, sS, size=n - 1) - sS ** 2 / 2)
            ia, pred = legs + svc, tb + sb
            m = (ia > 0) & (pred > 0)
            resid.append(np.log(ia[m] / pred[m]))
            grp.append(np.full(m.sum(), j))
        res[label] = re_decomp(np.concatenate(resid), np.concatenate(grp))
    return res


def main():
    root, out = sys.argv[1], sys.argv[2]
    info = dict(
        station_days_train=station_days(os.path.join(root, "data/raw/amazon_lmrrc_train/route_data.json")),
        station_days_eval=station_days(os.path.join(root, "data/raw/amazon_lmrrc_eval/eval_route_data.json")),
        konovalenko=konovalenko(os.path.join(root, "data/processed/konovalenko_stops.pkl")),
        simulated=simulated(os.path.join(root, "data/processed/train_routes.pkl")),
    )
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(info, open(out, "w"), indent=1)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
