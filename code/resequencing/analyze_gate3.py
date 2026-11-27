"""Gate 3 analysis. Implements the decision rules in analysis/gate3/preregistration.md.

Usage:
    python -I code/resequencing/analyze_gate3.py <gate3_A.csv> <gate3_B.csv> <out_dir> [--gate2c <gate2c_main.csv>]
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression

KS = (0.02, 0.05, 0.10)
LADDER = [("L1", "delay_plan_min"), ("L2", "L2_proj"), ("L3", "L3_prior"), ("L4f", "L4f"),
          ("L3m05", "L3_m05"), ("L4m05", "L4_m05"), ("L3m2", "L3_m2"), ("L4m2", "L4_m2")]
WINDOW_REGIMES = ["MIXB", "MIXB20", "MIXB50", "TW"]
REG_ORDER = ["ETA", "MIXA", "DAY", "MIXB", "MIXB20", "MIXB50", "TW"]


def topk(target, score, K):
    """Sum of target over the top-K entries of score (stable: ties keep route order)."""
    return float(target[np.argsort(-score, kind="stable")[:K]].sum())


def eps_star(D, rho):
    """min over nondecreasing g of max_r |rho_r - g(D_r)| = 1/2 max_{D_a <= D_b} (rho_a - rho_b)^+."""
    gb = pd.DataFrame({"D": D, "r": rho}).groupby("D")["r"].agg(["max", "min"]).sort_index()
    pm = np.maximum.accumulate(gb["max"].to_numpy())
    return 0.5 * float(np.max(pm - gb["min"].to_numpy()))


def mechanism(g, share=0.05):
    """Proposition 2 terms and losses in captured risk for one cell (g sorted by ridx)."""
    rho = g.rho_min.to_numpy()
    D = g.delay_plan_min.to_numpy()
    N = len(g)
    K = max(1, int(round(share * N)))
    v_opt = topk(rho, rho, K)
    loss_L1 = 1 - topk(rho, D, K) / v_opt
    eps = eps_star(D, rho)
    B_i = 2 * K * eps / v_opt
    iso = IsotonicRegression(increasing=True, out_of_bounds="clip").fit(D, rho).predict(D)
    H_iso = float(np.sum((rho - iso) ** 2) / np.sum((rho - rho.mean()) ** 2))
    out = dict(K=K, V_opt_rho=v_opt, loss_rho_L1=loss_L1, eps_star=eps, B_i=B_i, H_iso=H_iso,
               bound_i_holds=bool(loss_L1 <= B_i + 1e-9))
    if "tau_ii_min" in g and g.tau_ii_min.notna().all():
        L2 = g.L2_proj.to_numpy()
        loss_L2 = 1 - topk(rho, L2, K) / v_opt
        B_ii = 2 * K * float(g.tau_ii_min.max()) / v_opt
        out.update(loss_rho_L2=loss_L2, B_ii=B_ii, bound_ii_holds=bool(loss_L2 <= B_ii + 1e-9),
                   route_ineq_ii_holds=bool(np.all(np.abs(L2 - rho) <= g.tau_ii_min.to_numpy() + 1e-6)))
    return out


def cell_record(g):
    g = g.sort_values("ridx").reset_index(drop=True)
    d = g.delta_min.to_numpy()
    rho = g.rho_min.to_numpy()
    N = len(g)
    rec = dict(routes=N, window_share=g.window_share_route.mean(), late_routes=float((rho > 0.5).mean()),
               omega=g.b_late_window_stops.sum() / g.b_late_stops.sum() if g.b_late_stops.sum() > 0 else np.nan)
    for s in KS:
        K = max(1, int(round(s * N)))
        v4 = topk(d, rho, K)
        tag = int(s * 100)
        rec[f"L0_eff_{tag}"] = d.mean() * K / v4 if v4 > 0 else np.nan
        for lv, col in LADDER:
            if col in g and g[col].notna().all():
                rec[f"{lv}_eff_{tag}"] = topk(d, g[col].to_numpy(), K) / v4 if v4 > 0 else np.nan
        vo = topk(d, d, K)
        rec[f"L4_vs_oracle_{tag}"] = v4 / vo if vo > 0 else np.nan
    rec.update(mechanism(g))
    return rec


def within_day(g, share=0.10, fixed_k=None):
    """Within-station-day selection: each day selects max(1, round(share * N_d)) routes (or fixed_k)."""
    num = {lv: 0.0 for lv, _ in LADDER[:4]}
    den = 0.0
    for _, h in g.sort_values("ridx").groupby(["station", "date"]):
        K = fixed_k or max(1, int(round(share * len(h))))
        d = h.delta_min.to_numpy()
        den += topk(d, h.rho_min.to_numpy(), K)
        for lv, col in LADDER[:4]:
            num[lv] += topk(d, h[col].to_numpy(), K)
    return {f"{lv}_eff_day": (num[lv] / den if den > 0 else np.nan) for lv in num}


def cells_A(df):
    rows = []
    for (sU, reg, f, nz), g in df.groupby(["sU", "regime", "frac", "noise"]):
        rec = dict(sU=sU, regime=reg, frac=f, noise=nz)
        rec.update(cell_record(g))
        if sU == 0:
            rec.update(within_day(g))
            rec.update({k.replace("_day", "_day_k1"): v for k, v in within_day(g, fixed_k=1).items()})
        rows.append(rec)
    t = pd.DataFrame(rows)
    t["regime"] = pd.Categorical(t.regime, REG_ORDER, ordered=True)
    return t.sort_values(["sU", "regime", "frac", "noise"]).reset_index(drop=True)


def cells_B(df):
    rows = []
    for (nz, phi, reg), g in df.groupby(["noise", "phi", "regime"]):
        rec = dict(V_level=nz, phi=phi, regime=reg, sR=g.sR.iloc[0], sE=g.sE.iloc[0])
        rec.update(cell_record(g))
        rows.append(rec)
    t = pd.DataFrame(rows)
    t["regime"] = pd.Categorical(t.regime, REG_ORDER, ordered=True)
    return t.sort_values(["V_level", "regime", "phi"]).reset_index(drop=True)


def pick(t, **kw):
    m = np.ones(len(t), bool)
    for k, v in kw.items():
        m &= (t[k] == v).to_numpy() if not isinstance(v, float) else np.isclose(t[k].astype(float), v)
    x = t[m]
    assert len(x) == 1, (kw, len(x))
    return x.iloc[0]


def decide(A, B):
    L = []
    res = {}
    # H5a / H5b at sigma_U = 0.50, f = 0.35, K/N = 5%
    a = A[np.isclose(A.sU, 0.5) & np.isclose(A.frac, 0.35)]
    cand = a[a.regime.isin(WINDOW_REGIMES)]
    qual = cand[cand.omega >= 0.5]
    fails = qual[qual.L2_eff_5 < qual.L1_eff_5]
    L.append("H5a (sigma_U = 0.50, f = 0.35, K/N = 5%): cells with omega >= 0.5 and V(L2) >= V(L1)?")
    for _, c in cand.iterrows():
        q = c.omega >= 0.5
        L.append(f"   {c.regime:7s} {c.noise:4s} omega={c.omega:.3f} {'qualifies' if q else 'not qualifying'}; "
                 f"L1={c.L1_eff_5:.3f} L2={c.L2_eff_5:.3f}" + (f" -> {'pass' if c.L2_eff_5 >= c.L1_eff_5 else 'FAIL'}" if q else ""))
    nq, nf = len(qual), len(fails)
    h5a_all = nq > 0 and nf == 0
    h5a_third = nq > 0 and nf <= nq / 3
    L.append(f"H5a: {nq - nf}/{nq} qualifying cells pass")
    b_fail = a[a.L3_eff_5 < 0.80]
    h5b = len(a) == 12 and len(b_fail) == 0
    L.append("H5b (sigma_U = 0.50, f = 0.35): V(L3)/V(L4*) >= 0.80 in all 12 cells? "
             + ", ".join(f"{c.regime}/{c.noise}={c.L3_eff_5:.3f}" for _, c in a.iterrows()) + f" -> {'pass' if h5b else 'FAIL'}")
    if h5a_all and h5b:
        oc = "GO"
    elif h5a_all and not h5b:
        oc = "GO-NARROW (forecast)"
    elif h5a_third:
        oc = "GO-NARROW (projection)"
    else:
        oc = "REVISE"
    res.update(H5a_qualifying=int(nq), H5a_fails=int(nf), H5b=bool(h5b), outcome=oc)
    # H3: Arm B series MIXB50, TW x base, high
    met = 0
    for reg in ("MIXB50", "TW"):
        for lvl in ("base", "high"):
            s = B[(B.regime == reg) & (B.V_level == lvl)].sort_values("phi")
            ys = s.L1_eff_5.to_numpy()
            gain = ys[-1] - ys[0]
            rs = spearmanr(s.phi, ys).statistic
            ok = gain >= 0.10 and rs >= 0.90
            met += ok
            L.append(f"   H3 series {reg} {lvl}: L1 eff over phi {list(s.phi)} = " + ", ".join(f"{y:.3f}" for y in ys)
                     + f"; gain {gain:+.3f}, Spearman {rs:.2f} -> {'meets' if ok else 'does not meet'}")
    h3 = met >= 3
    L.append(f"H3: {met}/4 series meet both conditions (need >= 3) -> {'supported' if h3 else 'not supported'}")
    h3d_v = pick(B, V_level="base", phi=0.75, regime="MIXB").L1_eff_5
    h3d = h3d_v <= 0.80
    L.append(f"H3-D: phi = 0.75, base, MIXB, L1 eff = {h3d_v:.3f} (<= 0.80) -> {'met' if h3d else 'not met'}")
    # H4 on the 28 W cells
    w = A[A.sU == 0]
    s4 = spearmanr(w.B_i, w.loss_rho_L1).statistic
    s4b = spearmanr(w.B_i, 1 - w.L1_eff_5).statistic
    h4 = s4 >= 0.80
    L.append(f"H4: Spearman(B_i, captured-risk loss of L1) over {len(w)} W cells = {s4:.3f} (>= 0.80) -> {'met' if h4 else 'not met'}"
             f"; vs. 1 - L1 eff: {s4b:.3f}")
    L.append(f"   secondary: Spearman(H_iso, loss) = {spearmanr(w.H_iso, w.loss_rho_L1).statistic:.3f}; "
             f"Spearman(omega, loss) = {spearmanr(w.omega, w.loss_rho_L1).statistic:.3f}; "
             f"Spearman(B_ii, L2 loss) = {spearmanr(w.B_ii, w.loss_rho_L2).statistic:.3f}; "
             f"vacuous B_i >= 1 in {(w.B_i >= 1).mean():.0%} of cells, B_ii >= 1 in {(w.B_ii >= 1).mean():.0%}")
    # H6 within station-days
    h6v = pick(A, sU=0.0, regime="MIXB", frac=0.35, noise="base").L1_eff_day
    h6 = h6v <= 0.80
    L.append(f"H6: within-day L1 eff (MIXB, f = 0.35, base) = {h6v:.3f} (<= 0.80) -> {'met' if h6 else 'not met'}")
    # bound checks (code checks; must all hold)
    allc = pd.concat([A, B], ignore_index=True)
    chk = dict(bound_i=bool(allc.bound_i_holds.all()), bound_ii=bool(allc.bound_ii_holds.all()),
               route_ineq_ii=bool(allc.route_ineq_ii_holds.all()))
    L.append(f"Bound checks (all cells): Proposition 2(i) {chk['bound_i']}, 2(ii) {chk['bound_ii']}, route-level 2(ii) {chk['route_ineq_ii']}")
    L.append(f"Outcome: {oc}")
    res.update(H3_series_met=int(met), H3=bool(h3), H3D_value=float(h3d_v), H3D=bool(h3d), H4_spearman=float(s4),
               H4_spearman_delta=float(s4b), H4=bool(h4), H6_value=float(h6v), H6=bool(h6), checks=chk)
    return res, L


def gate2c_replication(path):
    df = pd.read_csv(path)
    rows = []
    for (reg, f, nz), g in df.groupby(["regime", "frac", "noise"]):
        g = g.sort_values("ridx").reset_index(drop=True)
        rec = dict(regime=reg, frac=f, noise=nz)
        rec.update(mechanism(g.drop(columns=[c for c in ("tau_ii_min",) if c in g])))
        K = max(1, int(round(0.05 * len(g))))
        v4 = topk(g.delta_min.to_numpy(), g.rho_min.to_numpy(), K)
        rec["L1_eff_5"] = topk(g.delta_min.to_numpy(), g.delay_plan_min.to_numpy(), K) / v4
        rec["omega"] = g.b_late_window_stops.sum() / g.b_late_stops.sum() if g.b_late_stops.sum() > 0 else np.nan
        rows.append(rec)
    t = pd.DataFrame(rows)
    s = spearmanr(t.B_i, t.loss_rho_L1).statistic
    return t, s


INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # categorical slots 1-4, fixed order (validated)
MARKS = ["o", "s", "^", "D"]


def _style(plt):
    plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})


def figures(A, B, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _style(plt)
    # Figure 1: ladder efficiency vs. prediction error sigma_U (f = 0.35, K/N = 5%)
    regs = ["DAY", "MIXB", "MIXB50", "TW"]
    series = [("L1", "Behind schedule (L1)"), ("L2", "Projection (L2)"), ("L3", "Prior forecast (L3)"),
              ("L4f", "Posterior forecast (L4)")]
    a = A[np.isclose(A.frac, 0.35)]
    fig, axes = plt.subplots(2, 4, figsize=(10.4, 4.9), sharey=True, sharex=True)
    for i, nz in enumerate(("base", "high")):
        for j, reg in enumerate(regs):
            ax = axes[i, j]
            g = a[(a.regime == reg) & (a.noise == nz)].sort_values("sU")
            for k, (lv, lab) in enumerate(series):
                ax.plot(g.sU, g[f"{lv}_eff_5"], lw=2, color=SLOTS[k], marker=MARKS[k], ms=4.5, mec="white", mew=0.6,
                        label=lab)
            om = g[np.isclose(g.sU, 0.5)].omega
            ttl = f"{reg}, {nz} noise" + (f"\nω = {om.iloc[0]:.2f} at σ_U = 0.5" if len(om) else "")
            ax.set_title(ttl, fontsize=8, color=INK, loc="left")
            ax.set_ylim(0, 1.05)
            ax.set_xticks([0, 0.25, 0.5, 1.0])
            ax.set_xticklabels(["0", "0.25", "0.50", "1.00"])
            ax.axvline(0.5, color=INK2, lw=0.6, ls=":")
            ax.grid(color=GRID, lw=0.6)
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="upper center", ncol=4, frameon=False, fontsize=8)
    fig.supxlabel("Prediction error σ_U (log-SD of true mean / carrier mean); dotted line = primary level 0.50",
                  fontsize=8.5, color=INK)
    fig.supylabel("Share of risk-triage value captured (K/N = 5%)", fontsize=8.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    plt.rcParams["pdf.fonttype"] = 42
    for ext, kw in (("png", {"dpi": 600}), ("pdf", {})):  # the manuscript uses the vector PDF
        fig.savefig(os.path.join(out_dir, f"gate3-prediction-error.{ext}"), **kw)
    plt.close(fig)
    # Figure 2 (main exhibit if H3 is supported): heat map of behind-schedule efficiency over
    # promise regime x route-wide share phi, at equal total variance
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("seq_blue", ["#eef4fc", "#9cc2ee", "#2a78d6", "#0f3f80"])
    phis = sorted(B.phi.unique())
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.0), sharey=True)
    for i, (lvl, phi0) in enumerate((("base", 0.0225 / 0.085), ("high", 0.0625 / 0.185))):
        ax = axes[i]
        M = np.array([[B[(B.regime == reg) & (B.V_level == lvl) & np.isclose(B.phi, p)].L1_eff_5.iloc[0]
                       for p in phis] for reg in regs])
        im = ax.imshow(M, cmap=cmap, vmin=0.2, vmax=1.0, aspect="auto")
        for r in range(M.shape[0]):
            for c in range(M.shape[1]):
                ax.text(c, r, f"{M[r, c]:.2f}", ha="center", va="center", fontsize=8,
                        color="white" if M[r, c] >= 0.70 else INK)
        c75 = phis.index(0.75)
        ax.add_patch(plt.Rectangle((c75 - 0.5, -0.5), 1, len(regs), fill=False, ec=INK, lw=1.4))
        ax.set_xticks(range(len(phis)))
        ax.set_xticklabels([f"{p:.2f}" for p in phis])
        ax.set_yticks(range(len(regs)))
        ax.set_yticklabels(regs)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xlabel("Route-wide share φ")
        ax.set_title(f"{lvl.capitalize()} total variance (assumed φ = {phi0:.2f})", fontsize=8.5, color=INK,
                     loc="left")
    cb = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cb.set_label("Behind schedule / risk triage (K/N = 5%)", fontsize=8, color=INK)
    cb.outline.set_visible(False)
    for ext, kw in (("png", {"dpi": 600}), ("pdf", {})):
        fig.savefig(os.path.join(out_dir, f"gate3-delay-structure.{ext}"), bbox_inches="tight", **kw)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv_a")
    ap.add_argument("csv_b")
    ap.add_argument("out_dir")
    ap.add_argument("--gate2c", default=None)
    ap.add_argument("--figures", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    A = cells_A(pd.read_csv(a.csv_a))
    B = cells_B(pd.read_csv(a.csv_b))
    A.to_csv(os.path.join(a.out_dir, "gate3_cells_A.csv"), index=False)
    B.to_csv(os.path.join(a.out_dir, "gate3_cells_B.csv"), index=False)
    cols = ["sU", "regime", "frac", "noise", "routes", "omega", "L0_eff_5", "L1_eff_5", "L2_eff_5", "L3_eff_5", "L4f_eff_5",
            "L4_vs_oracle_5", "B_i", "loss_rho_L1", "L1_eff_day"]
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(A[[c for c in cols if c in A]].round(3).to_string(index=False))
        print()
        print(B[["V_level", "phi", "regime", "omega", "L1_eff_5", "L2_eff_5", "L3_eff_5", "L4_vs_oracle_5", "B_i"]].round(3).to_string(index=False))
    res, lines = decide(A, B)
    print("\n" + "\n".join(lines))
    if a.gate2c:
        t2, s2 = gate2c_replication(a.gate2c)
        t2.to_csv(os.path.join(a.out_dir, "gate3_h4_replication_gate2c.csv"), index=False)
        res["H4_replication_gate2c_spearman"] = float(s2)
        res["H4_replication_bound_i_holds"] = bool(t2.bound_i_holds.all())
        print(f"\nH4 replication on Gate 2c cells: Spearman(B_i, loss) = {s2:.3f}; bound holds in all cells: {bool(t2.bound_i_holds.all())}")
    json.dump(res, open(os.path.join(a.out_dir, "gate3_decision.json"), "w"), indent=1)
    if a.figures:
        figures(A, B, a.out_dir)


if __name__ == "__main__":
    main()
