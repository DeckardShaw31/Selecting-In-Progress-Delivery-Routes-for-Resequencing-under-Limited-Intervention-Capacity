"""First-order route share of simulated inter-arrival residuals (manuscript Appendix B, equation route-share).

Compares the approximation
    share ~ phi * E[abar^2] sigma^2 / (E[a^2] sigma^2 + E[(1 - a)^2] sigma_S^2),
with a = tbar / (tbar + sbar) the travel share of a leg's planned inter-arrival time and abar its route mean,
against the simulated shares in analysis/gate3/calibration.json (gate3_calibrate.simulated: same 300 routes,
same seed). Exploratory check, added 2026-11-27; no triage outcome is computed.

Usage:
    python -I code/resequencing/route_share_check.py <project_root> <out_csv>
"""
from __future__ import annotations

import json
import os
import pickle
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate1 import NOISE  # noqa: E402
from gate3_calibrate import PHI_GRID  # noqa: E402


def leg_moments(train_pkl, n_routes=300, seed=3031):
    """Leg-weighted moments of the travel share a over the routes used by gate3_calibrate.simulated."""
    routes = pickle.load(open(train_pkl, "rb"))
    elig = [i for i, r in enumerate(routes) if r["n"] - 1 >= 10]
    pick = np.random.default_rng(seed).choice(elig, size=n_routes, replace=False)
    a_all, abar2, legs = [], [], []
    for i in pick:
        r = routes[i]
        n = r["n"]
        tb = np.array([float(r["tbar"][k - 1, k]) for k in range(1, n)])
        sb = r["sbar"][: n - 1].astype(float)  # service at the previous position, as in the simulation
        pred = tb + sb
        a = tb[pred > 0] / pred[pred > 0]
        a_all.append(a)
        abar2.append(a.mean() ** 2)
        legs.append(len(a))
    a_all, legs = np.concatenate(a_all), np.array(legs)
    return dict(mean_a=float(a_all.mean()), E_a2=float(np.mean(a_all ** 2)), E_1ma2=float(np.mean((1 - a_all) ** 2)),
                E_abar2=float(np.sum(legs * np.array(abar2)) / legs.sum()))


def main():
    root, out = sys.argv[1], sys.argv[2]
    mom = leg_moments(os.path.join(root, "data/processed/train_routes.pkl"))
    sim = json.load(open(os.path.join(root, "analysis/gate3/calibration.json")))["simulated"]

    def first_order(sR, sE, sS):
        return mom["E_abar2"] * sR ** 2 / (mom["E_a2"] * (sR ** 2 + sE ** 2) + mom["E_1ma2"] * sS ** 2)

    rows = []
    for lvl, (sR, sE, sS) in NOISE.items():
        rows.append(dict(setting=f"{lvl}_sU0.00", phi=sR ** 2 / (sR ** 2 + sE ** 2), first_order=first_order(sR, sE, sS),
                         simulated=sim[f"{lvl}_sU0.00"]["route_share"]))
    for lvl in ("base", "high"):
        sR0, sE0, sS0 = NOISE[lvl]
        V = sR0 ** 2 + sE0 ** 2
        for phi in PHI_GRID:
            rows.append(dict(setting=f"split_{lvl}_phi{phi:.2f}", phi=phi,
                             first_order=first_order(np.sqrt(phi * V), np.sqrt((1 - phi) * V), sS0),
                             simulated=sim[f"split_{lvl}_phi{phi:.2f}"]["route_share"]))
    d = pd.DataFrame(rows)
    d["ratio"] = d.first_order / d.simulated
    for k, v in mom.items():
        d[k] = v
    d.to_csv(out, index=False)
    print({k: round(v, 3) for k, v in mom.items()})
    print(d[["setting", "phi", "first_order", "simulated", "ratio"]].round(4).to_string(index=False))
    print(f"max relative deviation: {np.max(np.abs(d.ratio - 1)):.3f}")


if __name__ == "__main__":
    main()
