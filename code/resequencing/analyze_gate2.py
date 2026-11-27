"""Gate 2 analysis: information ladder (delay -> projection -> prior -> posterior risk) and the
positional-exposure mechanism. Implements the decision rule in analysis/gate2/preregistration.md.

Usage:
    python -I code/resequencing/analyze_gate2.py <gate2_results.csv> <out_dir> [--figures]
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

KS = [0.02, 0.05, 0.10]
LEVELS = [("L0", None), ("L1", "delay_plan_min"), ("L2", "L2_proj"), ("L3", "L3_prior"), ("L4", "rho_min")]
REGIME_ORDER = ["ETA", "MIXA", "MIXB", "TW"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SEQ = {"L1": "#86b6ef", "L2": "#5598e7", "L3": "#2a78d6", "L4": "#184f95"}  # one-hue ordinal ramp (blue 250-600)


def value(g, score_col, share, rng=None):
    N = len(g)
    K = max(1, int(round(share * N)))
    d = g["delta_min"].to_numpy()
    if score_col is None:
        return d.mean() * K  # expected value of random selection
    s = g[score_col].to_numpy()
    return d[np.argsort(-s, kind="stable")[:K]].sum()


def cell_table(df):
    rows = []
    for (reg, f, nz), g in df.groupby(["regime", "frac", "noise"]):
        b = g["b_late_stops"].to_numpy()
        cv_b = b.std() / b.mean() if b.mean() > 0 else np.nan
        bp = b[b > 0]
        rec = dict(regime=reg, frac=f, noise=nz, routes=len(g), late_routes=(g.rho_min > 0.5).mean(),
                   mean_rho=g.rho_min.mean(), R=g.delta_min.sum() / g.rho_min.sum() if g.rho_min.sum() > 0 else np.nan,
                   cv_b=cv_b, cv_b_pos=bp.std() / bp.mean() if len(bp) and bp.mean() > 0 else np.nan,
                   pred_eff=1.0 / np.sqrt(1.0 + cv_b ** 2) if np.isfinite(cv_b) else np.nan,
                   corr_D_rho=np.corrcoef(g.delay_plan_min, g.rho_min)[0, 1] if g.rho_min.std() > 0 else np.nan)
        for s in KS:
            v4 = value(g, "rho_min", s)
            vo = np.sort(g.delta_min.to_numpy())[::-1][:max(1, int(round(s * len(g))))].sum()
            for lv, col in LEVELS:
                rec[f"{lv}_eff_{int(s*100)}"] = value(g, col, s) / v4 if v4 > 0 else np.nan
            rec[f"L4_vs_oracle_{int(s*100)}"] = v4 / vo if vo > 0 else np.nan
            v1, v2, v3 = (value(g, c, s) for c in ("delay_plan_min", "L2_proj", "L3_prior"))
            gap = v4 - v1
            rec[f"gap_closed_L2_{int(s*100)}"] = (v2 - v1) / gap if gap > 1e-9 else np.nan
            rec[f"gap_closed_L3_{int(s*100)}"] = (v3 - v1) / gap if gap > 1e-9 else np.nan
        rows.append(rec)
    t = pd.DataFrame(rows)
    t["regime"] = pd.Categorical(t.regime, REGIME_ORDER, ordered=True)
    return t.sort_values(["regime", "frac", "noise"]).reset_index(drop=True)


def decision(t):
    p = t[(t.regime == "MIXB") & (t.frac == 0.35) & (t.noise == "base")].iloc[0]
    g2a = p.L1_eff_5 <= 0.70
    s1 = spearmanr(t.L1_eff_5, t.corr_D_rho).statistic
    s2 = spearmanr(t.corr_D_rho, t.pred_eff).statistic
    g2b = (s1 >= 0.70) and (s2 >= 0.70)
    outcome = "GO" if (g2a and g2b) else ("NARROW" if g2b else "REVISE")
    lines = [f"G2a: MIXB f=0.35 base K/N=5%: V_L1/V_L4 = {p.L1_eff_5:.3f} (threshold <= 0.70) -> {'pass' if g2a else 'fail'}",
             f"G2b(i): Spearman(V_L1/V_L4, corr(D,rho)) over {len(t)} cells = {s1:.3f} (>= 0.70)",
             f"G2b(ii): Spearman(corr(D,rho), 1/sqrt(1+CV(b)^2)) = {s2:.3f} (>= 0.70) -> G2b {'pass' if g2b else 'fail'}",
             f"Outcome: {outcome}"]
    return outcome, lines


def figures(t, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})
    # Figure 1: information ladder by regime (f=0.35, base, K/N=5%)
    sub = t[(t.frac == 0.35) & (t.noise == "base")].set_index("regime").loc[REGIME_ORDER]
    labels = {"L1": "L1 delay so far", "L2": "L2 + promise positions", "L3": "L3 + variability", "L4": "L4 + learning"}
    fig, ax = plt.subplots(figsize=(6.6, 3.0))
    x = np.arange(len(REGIME_ORDER))
    wbar = 0.19
    for k, lv in enumerate(["L1", "L2", "L3", "L4"]):
        ys = sub[f"{lv}_eff_5"].to_numpy()
        ax.bar(x + (k - 1.5) * (wbar + 0.02), ys, width=wbar, color=SEQ[lv], label=labels[lv])
        for xi, y in zip(x + (k - 1.5) * (wbar + 0.02), ys):
            ax.text(xi, y + 0.015, f"{y:.2f}", ha="center", fontsize=6.5, color=INK)
    ax.axhline(0.70, ls=(0, (4, 3)), lw=1, color=INK2)
    ax.text(-0.48, 0.715, "G2a threshold", color=INK2, fontsize=7.5, ha="left")
    ax.set_xticks(x, ["ETA\nwindows", "ETA + customer\nwindows (MIXA)", "Day-definite +\ncustomer windows\n(MIXB)", "Customer\nwindows only\n(TW)"])
    ax.set_ylabel("Lateness removed vs. L4 risk triage")
    ax.set_ylim(0, 1.12)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.legend(frameon=False, fontsize=7.5, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.16))
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "gate2-information-ladder.png"), dpi=200)
    plt.close(fig)
    # Figure 2: mechanism - delay-triage efficiency vs predicted efficiency 1/sqrt(1+CV(b)^2)
    fig, ax = plt.subplots(figsize=(3.6, 3.2))
    marks = {"ETA": "o", "MIXA": "s", "MIXB": "^", "TW": "D"}
    for reg in REGIME_ORDER:
        g = t[t.regime == reg]
        ax.scatter(g.pred_eff, g.L1_eff_5, s=28, marker=marks[reg], color="#2a78d6", edgecolor="white", linewidth=0.6, label=reg)
    ax.plot([0, 1], [0, 1], ls=(0, (4, 3)), lw=1, color=INK2)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Predicted efficiency 1/sqrt(1 + CV(b)^2)")
    ax.set_ylabel("Delay triage / risk triage (K/N = 5%)")
    ax.grid(color=GRID, lw=0.6)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right", title="promise regime", title_fontsize=7.5)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "gate2-exposure-mechanism.png"), dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("out_dir")
    ap.add_argument("--figures", action="store_true")
    a = ap.parse_args()
    df = pd.read_csv(a.csv)
    os.makedirs(a.out_dir, exist_ok=True)
    t = cell_table(df)
    t.to_csv(os.path.join(a.out_dir, "gate2_cells.csv"), index=False)
    cols = ["regime", "frac", "noise", "late_routes", "R", "cv_b", "pred_eff", "corr_D_rho",
            "L0_eff_5", "L1_eff_5", "L2_eff_5", "L3_eff_5", "gap_closed_L2_5", "gap_closed_L3_5", "L4_vs_oracle_5",
            "L1_eff_2", "L1_eff_10"]
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(t[cols].round(3).to_string(index=False))
    outcome, lines = decision(t)
    print("\n" + "\n".join(lines))
    if a.figures:
        figures(t, a.out_dir)


if __name__ == "__main__":
    main()
