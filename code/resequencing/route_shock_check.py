"""Channel-2 check of Proposition 7 on Gate 3 Arm B (exploratory; expectation in analysis/theory/proof_audit.md 3c).

In the DAY regime there is no waiting, so delay_plan_min is the sum of the deviations of the completed legs and
services. With lognormal legs t * exp(theta + eps - V/2), theta ~ N(0, phi V), eps ~ N(0, (1 - phi) V), and services
s * exp(zeta - sS^2/2), the exact variance is

    Var D = (e^{phi V} - 1)(Lam^2 - Q) + (e^V - 1) Q + (e^{sS^2} - 1) sum s^2,

with Lam = sum of completed planned leg times and Q = sum of their squares. The cross-route mean of D^2 in each DAY
cell should equal the mean predicted variance.

Usage:
    python -I code/resequencing/route_shock_check.py <train_routes.pkl> <sample_routes.csv> <gate3_B.csv> <out_dir>
"""
from __future__ import annotations

import json
import os
import pickle
import sys

import numpy as np
import pandas as pd

NOISE = {"base": (0.15, 0.25, 0.40), "high": (0.25, 0.35, 0.60)}
FRAC = 0.35


def route_terms(r):
    n = r["n"]
    m = int(round(FRAC * (n - 1)))
    if m < 1 or n - 1 - m < 2:
        return None
    t = np.array([float(r["tbar"][k - 1, k]) for k in range(1, m + 1)])
    s = np.asarray(r["sbar"], dtype=float)[1:m + 1]
    lam, q = t.sum(), (t ** 2).sum()
    return dict(Lam=lam, Q=q, n_eff=lam ** 2 / q if q > 0 else np.nan, legs=m, S2=(s ** 2).sum())


def main(pkl, sample, gate3b, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    routes = pickle.load(open(pkl, "rb"))
    ridx = pd.read_csv(sample).ridx.to_numpy()
    terms = pd.DataFrame({i: route_terms(routes[i]) for i in ridx}).T
    terms.index.name = "ridx"
    b = pd.read_csv(gate3b, usecols=["ridx", "noise", "phi", "regime", "delay_plan_min"])
    b = b[b.regime == "DAY"].join(terms, on="ridx")
    rows = []
    for (lvl, phi), g in b.groupby(["noise", "phi"]):
        sR0, sE0, sS = NOISE[lvl]
        V = sR0 ** 2 + sE0 ** 2
        pred = ((np.exp(phi * V) - 1) * (g.Lam ** 2 - g.Q) + (np.exp(V) - 1) * g.Q
                + (np.exp(sS ** 2) - 1) * g.S2) / 3600.0  # seconds^2 -> minutes^2
        d2 = g.delay_plan_min ** 2
        ratio = d2.mean() / pred.mean()
        se = d2.std() / np.sqrt(len(d2)) / pred.mean()
        travel_pred = (((np.exp(phi * V) - 1) * (g.Lam ** 2 - g.Q) + (np.exp(V) - 1) * g.Q) / 3600.0).mean()
        rows.append(dict(noise=lvl, phi=phi, routes=len(g), mean_D2_obs=d2.mean(), mean_D2_pred=pred.mean(),
                         ratio=ratio, se_ratio=se, z=(ratio - 1) / se, within_3se=bool(abs(ratio - 1) <= 3 * se),
                         sd_obs_min=float(np.sqrt(d2.mean())), sd_pred_min=float(np.sqrt(pred.mean())),
                         travel_share_pred=float(travel_pred / pred.mean())))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(out_dir, "route_shock_check.csv"), index=False)
    ne = terms.n_eff.astype(float)
    growth = {}
    for lvl in ("base", "high"):
        g = t[t.noise == lvl].set_index("phi")
        growth[lvl] = dict(obs=float(g.sd_obs_min[0.9] / g.sd_obs_min[0.1]),
                           pred=float(g.sd_pred_min[0.9] / g.sd_pred_min[0.1]))
    summ = dict(cells=len(t), all_within_3se=bool(t.within_3se.all()),
                n_eff_median=float(ne.median()), n_eff_q10=float(ne.quantile(0.1)), n_eff_q90=float(ne.quantile(0.9)),
                legs_median=float(terms.legs.astype(float).median()), sd_growth_phi_0p1_to_0p9=growth)
    json.dump(summ, open(os.path.join(out_dir, "route_shock_check.json"), "w"), indent=1)
    with pd.option_context("display.width", 200):
        print(t.round(4).to_string(index=False))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:5])
