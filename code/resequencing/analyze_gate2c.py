"""Gate 2c analysis (confirmatory, training split). Implements analysis/gate2c/preregistration.md.

Usage:
    python -I code/resequencing/analyze_gate2c.py <gate2c_results.csv> <out_dir> [--figures]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_gate2 import value  # noqa: E402

DOSE = ["DAY", "MIXB", "MIXB20", "MIXB50"]
ALL = ["ETA", "MIXA", "DAY", "MIXB", "MIXB20", "MIXB50", "TW"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SERIES = {("0.35", "base"): "#2a78d6", ("0.5", "base"): "#eb6834", ("0.35", "high"): "#1baf7a", ("0.5", "high"): "#e87ba4"}


def cells(df):
    rows = []
    for (reg, f, nz), g in df.groupby(["regime", "frac", "noise"]):
        v4 = {s: value(g, "rho_min", s) for s in (0.02, 0.05, 0.10)}
        rec = dict(regime=reg, frac=f, noise=nz, routes=len(g),
                   window_share=g.window_share_route.mean(), late_routes=(g.rho_min > 0.5).mean(),
                   corr_b_D=np.corrcoef(g.b_late_stops, g.delay_plan_min)[0, 1],
                   window_exposure_share=g.b_late_window_stops.sum() / g.b_late_stops.sum() if g.b_late_stops.sum() > 0 else np.nan)
        for s in (0.02, 0.05, 0.10):
            for lv, col in [("L0", None), ("L1", "delay_plan_min"), ("L2", "L2_proj"), ("L3", "L3_prior")]:
                rec[f"{lv}_eff_{int(s*100)}"] = value(g, col, s) / v4[s] if v4[s] > 0 else np.nan
            vo = np.sort(g.delta_min.to_numpy())[::-1][:max(1, int(round(s * len(g))))].sum()
            rec[f"L4_vs_oracle_{int(s*100)}"] = v4[s] / vo if vo > 0 else np.nan
        rows.append(rec)
    t = pd.DataFrame(rows)
    t["regime"] = pd.Categorical(t.regime, ALL, ordered=True)
    return t.sort_values(["regime", "frac", "noise"]).reset_index(drop=True)


def decide(t):
    out = []
    c1v = t[(t.regime == "MIXB") & (t.frac == 0.35) & (t.noise == "base")].L1_eff_5.iloc[0]
    c1 = c1v <= 0.70
    out.append(f"C1: MIXB f=0.35 base K/N=5%: V_L1/V_L4 = {c1v:.3f} (<= 0.70) -> {'pass' if c1 else 'fail'}")
    mono = 0
    for f in (0.35, 0.5):
        for nz in ("base", "high"):
            ys = [t[(t.regime == r) & (t.frac == f) & (t.noise == nz)].L1_eff_5.iloc[0] for r in DOSE]
            ok = all(ys[k + 1] <= ys[k] for k in range(3))
            mono += ok
            out.append(f"   dose series f={f} {nz}: " + " -> ".join(f"{y:.3f}" for y in ys) + f"  {'monotone' if ok else 'NOT monotone'}")
    d0 = t[(t.regime == "DAY") & (t.frac == 0.35) & (t.noise == "base")].L1_eff_5.iloc[0]
    d50 = t[(t.regime == "MIXB50") & (t.frac == 0.35) & (t.noise == "base")].L1_eff_5.iloc[0]
    c2 = (mono >= 3) and (d0 - d50 >= 0.20)
    out.append(f"C2: monotone in {mono}/4 combos (need >= 3); DAY - MIXB50 drop = {d0 - d50:.3f} (need >= 0.20) -> {'pass' if c2 else 'fail'}")
    s3 = spearmanr(t.L1_eff_5, t.corr_b_D).statistic
    c3 = s3 >= 0.70
    out.append(f"C3: Spearman(V_L1/V_L4, corr(b,D)) over {len(t)} cells = {s3:.3f} (>= 0.70) -> {'pass' if c3 else 'fail'}")
    if c1 and c2 and c3:
        oc = "GO"
    elif c2 and c1 and not c3:
        oc = "GO-NARROW (mechanism statistic weak)"
    elif c2 and not c1:
        oc = "GO-NARROW (gap smaller on fresh routes; regime map is the main result)"
    else:
        oc = "STOP option 1 as a TR-E paper (dose-response not confirmed)"
    out.append(f"Outcome: {oc}")
    return oc, out


def figures(t, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.2))
    ax = axes[0]
    for (f, nz), col in SERIES.items():
        g = t[(t.frac == float(f)) & (t.noise == nz) & (t.regime.isin(DOSE))].set_index("regime").loc[DOSE]
        xs = g.window_share.to_numpy() * 100
        ax.plot(xs, g.L1_eff_5, lw=2, color=col, marker="o", ms=4, label=f"check-in {f}, {nz} noise")
    ax.set_xlabel("Stops with a customer-window promise (%)")
    ax.set_ylabel("Delay triage / risk triage (K/N = 5%)")
    ax.set_ylim(0, 1.05)
    ax.grid(color=GRID, lw=0.6)
    ax.legend(frameon=False, fontsize=7, loc="lower left")
    ax.set_title("Dose-response (training split)", fontsize=9, color=INK, loc="left")
    ax = axes[1]
    marks = {"ETA": "o", "MIXA": "s", "DAY": "P", "MIXB": "^", "MIXB20": "v", "MIXB50": "X", "TW": "D"}
    for reg in ALL:
        g = t[t.regime == reg]
        ax.scatter(g.corr_b_D, g.L1_eff_5, s=26, marker=marks[reg], color="#2a78d6", edgecolor="white", linewidth=0.6, label=reg)
    ax.set_xlabel("Exposure endogeneity corr(b, D)")
    ax.set_ylim(0, 1.05)
    ax.grid(color=GRID, lw=0.6)
    ax.legend(frameon=False, fontsize=6.5, loc="lower right", ncol=2)
    ax.set_title("Pre-specified mechanism statistic (failed)", fontsize=9, color=INK, loc="left")
    ax = axes[2]
    for reg in ALL:
        g = t[t.regime == reg]
        ax.scatter(g.window_exposure_share, g.L1_eff_5, s=26, marker=marks[reg], color="#2a78d6", edgecolor="white", linewidth=0.6)
    ax.set_xlabel("At-risk promises that are windows (share)")
    ax.set_ylim(0, 1.05)
    ax.set_xlim(-0.03, 1.03)
    ax.grid(color=GRID, lw=0.6)
    ax.set_title("Secondary measure (no threshold)", fontsize=9, color=INK, loc="left")
    fig.tight_layout()
    plt.rcParams["pdf.fonttype"] = 42
    for ext, kw in (("png", {"dpi": 600}), ("pdf", {})):  # the manuscript uses the vector PDF
        fig.savefig(os.path.join(out_dir, f"gate2c-dose-response-mechanism.{ext}"), **kw)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("out_dir")
    ap.add_argument("--figures", action="store_true")
    a = ap.parse_args()
    df = pd.read_csv(a.csv)
    os.makedirs(a.out_dir, exist_ok=True)
    t = cells(df)
    t.to_csv(os.path.join(a.out_dir, "gate2c_cells.csv"), index=False)
    cols = ["regime", "frac", "noise", "routes", "window_share", "late_routes", "corr_b_D", "window_exposure_share",
            "L0_eff_5", "L1_eff_5", "L2_eff_5", "L3_eff_5", "L4_vs_oracle_5", "L1_eff_2", "L1_eff_10"]
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(t[cols].round(3).to_string(index=False))
    _, lines = decide(t)
    print("\n" + "\n".join(lines))
    if a.figures:
        figures(t, a.out_dir)


if __name__ == "__main__":
    main()
