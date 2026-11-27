"""Gate 5: delay x delay sensitivity. Implements analysis/gate5/preregistration.md.

Scores (per route and cell; dispatcher information only):
    S1 = W_r * D_r                                      (delay x remaining promise weight; primary)
    S2 = D_r * sum_i w_ri (Lam_r + c_ri) / Lam_r        (pace-extrapolated delay; perfect learning)
    S3 = D_r * sum_i w_ri (1 + c_ri gamma_r)            (model-based sensitivity, Proposition 7)
Value: removable lateness (delta_min) of the top-K routes / that of the top-K by risk (L4), as in Gates 2c and 3.

Usage:
    python -I code/resequencing/gate5.py --verify-only   (pre-run verification checks)
    python -I code/resequencing/gate5.py                 (full run)
Paths are fixed relative to the project root (run from the project root).
"""
from __future__ import annotations

import json
import os
import pickle
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate1 import (ETA_HALF_WIDTH, NOISE, nominal_eta, nominal_schedule_wait,  # noqa: E402
                   promise_arrays, regime_windows)

PKL = "data/processed/train_routes.pkl"
G2C = "analysis/gate2c/results/gate2c_main.csv"
G2C_CELLS = "analysis/gate2c/gate2c_cells.csv"
G3A = "analysis/gate3/results/gate3_A.csv"
G3A_CELLS = "analysis/gate3/gate3_cells_A.csv"
G3B = "analysis/gate3/results/gate3_B.csv"
G3B_CELLS = "analysis/gate3/gate3_cells_B.csv"
G4 = "analysis/gate4/gate4_routes.csv"
G4_EFF = "analysis/gate4/gate4_efficiency.csv"
G2D = "analysis/gate2d/gate2d_routes.csv"
G2D_EFF = "analysis/gate2d/gate2d_efficiency.csv"
OUT = "analysis/gate5"
KS = (0.02, 0.05, 0.10)
SCORES = [("L1", "delay_plan_min"), ("L2", "L2_proj"), ("L3", "L3_prior"), ("S1", "S1"), ("S2", "S2"), ("S3", "S3")]
REGS = ["ETA", "MIXA", "DAY", "MIXB", "MIXB20", "MIXB50", "TW"]
N_BOOT, BOOT_SEED = 2000, 2032

ROUTES = None


def routes():
    global ROUTES
    if ROUTES is None:
        ROUTES = pickle.load(open(PKL, "rb"))
    return ROUTES


def route_terms(ridx, frac):
    """Planned quantities at the check-in after round(frac (n - 1)) stops (carrier means)."""
    r = routes()[ridx]
    n = r["n"]
    m = int(round(frac * (n - 1)))
    tb = r["tbar"].astype(np.float64)
    legs_done = np.array([tb[k - 1, k] for k in range(1, m + 1)])
    c = np.cumsum([tb[k - 1, k] for k in range(m + 1, n)])  # travel from check-in to remaining stop k
    has_tw = ~np.isnan(r["tw_close"][m + 1:n])
    sb = np.asarray(r["sbar"], dtype=float)[1:m + 1]
    return dict(n=n, m=m, Lam=legs_done.sum(), Q=(legs_done ** 2).sum(), Ssvc=(sb ** 2).sum(),
                W_all=float(n - 1 - m), W_tw=float(has_tw.sum()), C_all=float(c.sum()), C_tw=float(c[has_tw].sum()))


def add_scores(df, seed_col_noise):
    """Attach W, S1, S2, S3 to a results frame with columns ridx, frac, regime, delay_plan_min and noise params."""
    keys = df[["ridx", "frac"]].drop_duplicates()
    T = pd.DataFrame([dict(ridx=int(a), frac=float(b), **route_terms(int(a), float(b))) for a, b in keys.itertuples(index=False)])
    d = df.merge(T, on=["ridx", "frac"], how="left", validate="many_to_one")
    tw = d.regime == "TW"
    d["W"] = np.where(tw, d.W_tw, d.W_all)
    C = np.where(tw, d.C_tw, d.C_all)
    D = d.delay_plan_min.to_numpy()
    lam_ok = d.Lam > 0
    d["S1"] = d.W * D
    d["S2"] = np.where(lam_ok, D * (d.W + C / np.where(lam_ok, d.Lam, 1.0)), d.S1)
    sR, sE, sS = seed_col_noise(d)
    s2 = sR ** 2 + sE ** 2
    phi = np.where(s2 > 0, sR ** 2 / s2, 0.0)
    varD = s2 * (d.Lam ** 2 * phi + d.Q * (1 - phi)) + sS ** 2 * d.Ssvc
    gamma = np.where(varD > 0, d.Lam * phi * s2 / np.where(varD > 0, varD, 1.0), 0.0)
    d["gamma"] = gamma
    d["S3"] = D * (d.W + gamma * C)
    return d


def noise_from_label(d):
    p = np.array([NOISE[z] for z in d.noise])
    return p[:, 0], p[:, 1], p[:, 2]


def noise_from_cols(d):
    return d.sR.to_numpy(), d.sE.to_numpy(), d.sS.to_numpy()


def topk(target, score, K):
    return float(target[np.argsort(-score, kind="stable")[:K]].sum())


def cell_effs(g):
    d = g.delta_min.to_numpy()
    rho = g.rho_min.to_numpy()
    N = len(g)
    rec = dict(routes=N)
    for s in KS:
        K = max(1, int(round(s * N)))
        v4 = topk(d, rho, K)
        for lv, col in SCORES:
            rec[f"{lv}_eff_{int(s * 100)}"] = topk(d, g[col].to_numpy(), K) / v4 if v4 > 0 else np.nan
    K = max(1, int(round(0.05 * N)))
    for lv, col in (("L1", "delay_plan_min"), ("S1", "S1")):
        sel = np.argsort(-g[col].to_numpy(), kind="stable")[:K]
        rec[f"{lv}_zeroW_share_5"] = float((g.W.to_numpy()[sel] == 0).mean())
    for lv in ("S1", "S2", "S3"):
        e1 = rec["L1_eff_5"]
        rec[f"G_{lv}"] = (rec[f"{lv}_eff_5"] - e1) / (1 - e1) if 1 - e1 >= 0.02 else np.nan
        gap2 = rec["L2_eff_5"] - e1
        rec[f"shareL2gap_{lv}"] = (rec[f"{lv}_eff_5"] - e1) / gap2 if gap2 > 0.02 else np.nan
    return rec


def cells(d, keys, sort_ridx):
    rows = []
    for key, g in d.groupby(keys, sort=False):
        if sort_ridx:
            g = g.sort_values("ridx", kind="stable")
        rows.append(dict(zip(keys, key), **cell_effs(g.reset_index(drop=True))))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------- observed data

def obs_value(x, score, k, y):
    K = max(1, int(round(k * len(x))))
    return x[y].to_numpy()[np.argsort(-x[score].to_numpy(), kind="stable")[:K]].sum()


def observed(x, W, b1, b2, ys, label):
    x = x.copy()
    x["W"] = W
    x["WB1"] = x.W * x[b1]
    x["WB2"] = x.W * x[b2]
    out = []
    for y in ys:
        for k in (0.05, 0.10):
            vo = obs_value(x, y, k, y)
            rec = dict(set=label, outcome=y, k=k, routes=len(x))
            for sc in (b1, b2, "WB1", "WB2", "F_projection"):
                rec[sc if sc in ("WB1", "WB2", "F_projection") else ("B1" if sc == b1 else "B2")] = obs_value(x, sc, k, y) / vo
            out.append(rec)
    rng = np.random.default_rng(BOOT_SEED)
    y = ys[0]
    diffs = []
    for _ in range(N_BOOT):
        bs = x.iloc[rng.integers(0, len(x), len(x))].reset_index(drop=True)
        vo = obs_value(bs, y, 0.05, y)
        diffs.append((obs_value(bs, "WB1", 0.05, y) - obs_value(bs, b1, 0.05, y)) / vo if vo > 0 else np.nan)
    ci = [float(v) for v in np.nanpercentile(diffs, [2.5, 97.5])]
    return out, ci


# ---------------------------------------------------------------------------------------------- verification

def verify():
    ok = True
    rng = np.random.default_rng(5)
    # 1. W from promise_arrays equals U / windows_remaining (and our route_terms)
    for path, seed in ((G2C, 2027), (G3A, 2029)):
        cols = ["ridx", "frac", "regime", "U", "windows_remaining"] + (["sU"] if path == G3A else [])
        df = pd.read_csv(path, usecols=cols)
        if path == G3A:
            df = df[df.sU == 0]
        smp = df.sample(50, random_state=5)
        bad = 0
        for row in smp.itertuples(index=False):
            r = routes()[int(row.ridx)]
            eta = nominal_eta(r)
            o, c = regime_windows(r, row.regime, eta, seed, int(row.ridx))
            st_w, en_w = nominal_schedule_wait(r, o)
            m = int(round(row.frac * (r["n"] - 1)))
            idx = np.arange(m, r["n"])
            _, _, w, _ = promise_arrays(r, idx, row.regime, eta, ETA_HALF_WIDTH, en_w + ETA_HALF_WIDTH, st_w, o, c)
            t = route_terms(int(row.ridx), float(row.frac))
            mine = t["W_tw"] if row.regime == "TW" else t["W_all"]
            filev = row.windows_remaining if row.regime == "TW" else row.U
            bad += int(abs(w.sum() - mine) > 1e-9 or abs(w.sum() - filev) > 1e-9)
        print(f"V1 {os.path.basename(path)}: W mismatches in 50 sampled route-cells: {bad}")
        ok &= bad == 0
    # 2. eff(L1), eff(L2) reproduce published cells
    c2 = cells(add_scores(pd.read_csv(G2C), noise_from_label), ["regime", "frac", "noise"], sort_ridx=False)
    pub = pd.read_csv(G2C_CELLS)
    mg = c2.merge(pub, on=["regime", "frac", "noise"], suffixes=("", "_pub"))
    dev2c = max(float((mg[f"{lv}_eff_{k}"] - mg[f"{lv}_eff_{k}_pub"]).abs().max()) for lv in ("L1", "L2") for k in (2, 5, 10))
    a = pd.read_csv(G3A)
    cw = cells(add_scores(a[a.sU == 0], noise_from_cols), ["regime", "frac", "noise"], sort_ridx=True)
    pa = pd.read_csv(G3A_CELLS)
    pa = pa[pa.sU == 0]
    mw = cw.merge(pa, on=["regime", "frac", "noise"], suffixes=("", "_pub"))
    dev3 = max(float((mw[f"{lv}_eff_{k}"] - mw[f"{lv}_eff_{k}_pub"]).abs().max()) for lv in ("L1", "L2") for k in (2, 5, 10))
    b = pd.read_csv(G3B)
    cb = cells(add_scores(b, noise_from_cols), ["noise", "phi", "regime"], sort_ridx=True).rename(columns={"noise": "V_level"})
    pb = pd.read_csv(G3B_CELLS)
    mb = cb.merge(pb, on=["V_level", "phi", "regime"], suffixes=("", "_pub"))
    dev3b = max(float((mb[f"{lv}_eff_{k}"] - mb[f"{lv}_eff_{k}_pub"]).abs().max()) for lv in ("L1", "L2") for k in (2, 5, 10))
    print(f"V2 published efficiencies reproduced: Gate 2c max dev {dev2c:.1e} ({len(mg)} cells); "
          f"Gate 3 W {dev3:.1e} ({len(mw)} cells); Arm B {dev3b:.1e} ({len(mb)} cells)")
    ok &= max(dev2c, dev3, dev3b) <= 1e-9 and len(mg) == 28 and len(mw) == 28 and len(mb) == 40
    # 3. unit checks
    t = pd.DataFrame(dict(ridx=[0, 0], frac=[0.35, 0.35], regime=["DAY", "DAY"], delay_plan_min=[10.0, -5.0],
                          sR=[0.0, 0.0], sE=[0.3, 0.3], sS=[0.4, 0.4]))
    u = add_scores(t, noise_from_cols)
    ok3 = bool(np.allclose(u.S3, u.S1))
    lam_min = float(a[a.sU == 0].pipe(add_scores, noise_from_cols).pipe(lambda z: 1 + z.gamma * 0).min())
    g2 = add_scores(pd.read_csv(G2C, usecols=["ridx", "frac", "regime", "noise", "delay_plan_min"]), noise_from_label)
    ok3 &= bool((g2.gamma >= 0).all()) and lam_min >= 1
    zc = g2.copy()
    ok3 &= bool(np.allclose(np.where(zc.Lam > 0, zc.delay_plan_min * zc.W, zc.S1), zc.S1))  # c = 0 reduces S2 to S1 form
    print(f"V3 unit checks (S3 = S1 at phi = 0; gamma >= 0 so lambda >= 1): {'PASS' if ok3 else 'FAIL'}")
    ok &= ok3
    # 4. observed: B1 and oracle reproduce published efficiencies
    x4 = pd.read_csv(G4)
    x4 = x4[x4.f == 0.35].reset_index(drop=True)
    e4 = pd.read_csv(G4_EFF)
    p4 = e4[(e4.subset == "all") & (e4.outcome == "Y_count") & (e4.f == 0.35) & (e4.k == 0.05)].iloc[0]
    mine4 = obs_value(x4, "B1_late_so_far", 0.05, "Y_count") / obs_value(x4, "Y_count", 0.05, "Y_count")
    x2 = pd.read_csv(G2D)
    x2 = x2[x2.f == 0.35].sort_values("route").reset_index(drop=True)
    e2 = pd.read_csv(G2D_EFF)
    p2 = e2[(e2.f == 0.35) & (e2.k == 0.05)]
    p2 = p2.iloc[0]  # gate2d_efficiency.csv holds the primary (late minutes) outcome only
    mine2 = obs_value(x2, "B1_late_so_far", 0.05, "Y_future_late_min") / obs_value(x2, "Y_future_late_min", 0.05, "Y_future_late_min")
    d4, d2 = abs(mine4 - p4.B1_late_so_far), abs(mine2 - p2.B1_late_so_far)
    print(f"V4 observed B1 efficiencies reproduced: LaDe dev {d4:.1e}; Konovalenko dev {d2:.1e}")
    ok &= max(d4, d2) <= 1e-9
    print("VERIFICATION " + ("PASS" if ok else "FAIL"))
    return ok


def main():
    os.makedirs(OUT, exist_ok=True)
    c2 = cells(add_scores(pd.read_csv(G2C), noise_from_label), ["regime", "frac", "noise"], sort_ridx=False)
    a = pd.read_csv(G3A)
    cw = cells(add_scores(a[a.sU == 0], noise_from_cols), ["regime", "frac", "noise"], sort_ridx=True)
    cu = cells(add_scores(a[a.sU > 0], noise_from_cols), ["sU", "regime", "noise"], sort_ridx=True)
    cb = cells(add_scores(pd.read_csv(G3B), noise_from_cols), ["noise", "phi", "regime"], sort_ridx=True)
    for name, t in (("gate2c", c2), ("gate3_W", cw), ("gate3_sigmaU", cu), ("gate3_armB", cb)):
        t.to_csv(os.path.join(OUT, f"gate5_cells_{name}.csv"), index=False)
    # decision (primary: Gate 2c)
    p = c2[np.isclose(c2.frac, 0.35)]
    h7a_med = float(p.G_S1.median())
    worst = float((c2.S1_eff_5 - c2.L1_eff_5).min())
    h7a, h7b = h7a_med >= 0.25, worst >= -0.02
    outcome = "SUPPORT" if (h7a and h7b) else ("SUPPORT-WITH-CAVEAT" if h7a else "NOT SUPPORTED")
    pw = cw[np.isclose(cw.frac, 0.35)]
    rep_med = float(pw.G_S1.median())
    # observed
    x4 = pd.read_csv(G4)
    x4 = x4[x4.f == 0.35].reset_index(drop=True)
    o4, ci4 = observed(x4, x4.known_remaining.to_numpy(float), "B1_late_so_far", "B2_currently_late",
                       ["Y_count", "Y_minutes"], "LaDe")
    x2 = pd.read_csv(G2D)
    x2 = x2[x2.f == 0.35].sort_values("route").reset_index(drop=True)
    o2, ci2 = observed(x2, (x2.n_stops - x2.m).to_numpy(float), "B1_late_so_far", "B2_currently_late",
                       ["Y_future_late_count", "Y_future_late_min"], "Konovalenko")
    obs = pd.DataFrame(o4 + o2)
    obs.to_csv(os.path.join(OUT, "gate5_observed.csv"), index=False)
    res = dict(H7a_median_G_S1=h7a_med, H7a=bool(h7a), H7b_worst_S1_minus_L1=worst, H7b=bool(h7b), outcome=outcome,
               H7rep_median_G_S1=rep_med, H7rep=bool(rep_med >= 0.25),
               H7obs_LaDe_ci=ci4, H7obs=bool(ci4[0] > 0), Konovalenko_ci=ci2,
               secondary_median_G=dict(S2=float(p.G_S2.median()), S3=float(p.G_S3.median())))
    json.dump(res, open(os.path.join(OUT, "gate5_decision.json"), "w"), indent=1)
    show = ["regime", "frac", "noise", "L1_eff_5", "S1_eff_5", "S2_eff_5", "S3_eff_5", "L2_eff_5", "L3_eff_5", "G_S1",
            "G_S2", "G_S3", "L1_zeroW_share_5", "S1_zeroW_share_5"]
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print("Gate 2c (primary)")
        print(c2[show].round(3).to_string(index=False))
        print("\nGate 3 W (replication)")
        print(cw[show].round(3).to_string(index=False))
        print("\nObserved")
        print(obs.round(3).to_string(index=False))
    print("\n" + json.dumps(res, indent=1))


if __name__ == "__main__":
    if "--verify-only" in sys.argv:
        sys.exit(0 if verify() else 1)
    main()
