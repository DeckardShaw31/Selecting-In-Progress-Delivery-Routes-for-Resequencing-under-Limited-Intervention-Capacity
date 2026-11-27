"""Summarize Gate 1 results and draw the two Gate 1 figures.

Usage:
    python -I code/resequencing/analyze_gate1.py <results.csv> <out_dir> [--tag NAME] [--figures]
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

KS = [0.05, 0.10, 0.20]
SERIES = {"low": "#2a78d6", "base": "#eb6834", "high": "#1baf7a"}  # validated slots 1-3 (light)
MARKERS = {"low": "o", "base": "s", "high": "^"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def topk_stats(g, share):
    N = len(g)
    K = max(1, int(round(share * N)))
    d = g["delta_min"].to_numpy()
    rho = g["rho_min"].to_numpy()
    rec = np.argsort(-d, kind="stable")[:K]
    risk = np.argsort(-rho, kind="stable")[:K]
    v_rec, v_risk = d[rec].sum(), d[risk].sum()
    return v_risk / v_rec if v_rec > 0 else np.nan, len(set(rec) & set(risk)) / K


def summarize(df):
    rows = []
    for key, g in df.groupby(["regime", "frac", "noise"]):
        pos = g[g.rho_min > 0]
        R = pos.delta_min.sum() / pos.rho_min.sum() if len(pos) else np.nan
        ratio = (pos.delta_min / pos.rho_min).clip(lower=0) if len(pos) else pd.Series(dtype=float)
        rec = dict(regime=key[0], frac=key[1], noise=key[2], routes=len(g),
                   late_routes_share=(g.rho_min > 0.5).mean(),
                   mean_rho_min=g.rho_min.mean(), mean_delta_min=g.delta_min.mean(),
                   R_recoverable_share=R,
                   median_delta_over_rho=ratio[pos.rho_min > 0.5].median() if (pos.rho_min > 0.5).any() else np.nan,
                   share_routes_delta_ge_25pct=((pos.rho_min > 0.5) & (ratio >= 0.25)).sum() / max((pos.rho_min > 0.5).sum(), 1),
                   negative_delta_share=(g.delta_min < -1e-9).mean(),
                   spearman_rho_delta=spearmanr(g.rho_min, g.delta_min).statistic if g.delta_min.std() > 0 else np.nan)
        for s in KS:
            r, ov = topk_stats(g, s)
            rec[f"Vrisk_over_Vrec_{int(s*100)}"] = r
            rec[f"overlap_{int(s*100)}"] = ov
        rows.append(rec)
    return pd.DataFrame(rows)


def gate_decision(summ):
    b = summ[(summ.frac == 0.35) & (summ.noise == "base")]
    lines = []
    go = False
    for _, r in b.iterrows():
        g1a = r.R_recoverable_share >= 0.10
        g1b = r.Vrisk_over_Vrec_10 <= 0.85
        lines.append(f"- {r.regime}: R = {r.R_recoverable_share:.3f} (G1a {'pass' if g1a else 'fail'}); "
                     f"V_risk/V_rec at 10% = {r.Vrisk_over_Vrec_10:.3f} (G1b {'pass' if g1b else 'fail'})")
        go = go or (g1a and g1b)
    any_a = any(b.R_recoverable_share >= 0.10)
    outcome = "GO" if go else ("REVISE (risk ranking near-optimal)" if any_a else "REVISE (add a lever)")
    return outcome, lines


def figures(df, out_dir, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})
    regimes = [r for r in ["ETA", "TW"] if r in df.regime.unique()]
    # Figure 1: rho vs delta (base noise, f = 0.35)
    fig, axes = plt.subplots(1, len(regimes), figsize=(3.4 * len(regimes), 3.0), squeeze=False)
    unit = "late stops" if "objective" in df and (df.objective == "count").all() else "min"
    for ax, reg in zip(axes[0], regimes):
        g = df[(df.regime == reg) & (df.noise == "base") & (df.frac == 0.35)]
        g = g[g.rho_min > (0.05 if unit == "late stops" else 0.5)]
        lo = 0.01 if unit == "late stops" else 0.1
        mx = max(g.rho_min.max() if len(g) else 1, 1) * 1.3
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.plot([lo, mx], [lo, mx], ls=(0, (4, 3)), lw=1, color=INK2)
        xt = lo * (mx / lo) ** 0.08
        ax.text(xt, xt * (mx / lo) ** 0.45, "dashed: fully recoverable\nfloor: zero or negative", color=INK2, fontsize=7.5)
        ax.scatter(g.rho_min, g.delta_min.clip(lower=lo), s=12, color=SERIES["low"], edgecolor="white", linewidth=0.5)
        ax.set_xlim(lo, mx)
        ax.set_ylim(lo * 0.8, mx)
        ax.set_title(f"{reg} promises ({len(g)} late routes)", fontsize=9, color=INK, loc="left")
        ax.set_xlabel(f"Expected lateness, current sequence ({unit})")
        ax.grid(color=GRID, lw=0.6)
    axes[0][0].set_ylabel(f"Recoverable lateness ({unit}, log)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"gate1-rho-vs-delta{tag}.png"), dpi=200)
    plt.close(fig)
    # Figure 2: capacity curve V_risk / V_rec vs K/N (f = 0.35)
    shares = [0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50]
    fig, axes = plt.subplots(1, len(regimes), figsize=(3.4 * len(regimes), 3.0), squeeze=False, sharey=True)
    for ax, reg in zip(axes[0], regimes):
        ends = []
        for nz in ["low", "base", "high"]:
            g = df[(df.regime == reg) & (df.noise == nz) & (df.frac == 0.35)]
            if not len(g):
                continue
            ys = [topk_stats(g, s)[0] for s in shares]
            ax.plot([s * 100 for s in shares], ys, lw=2, color=SERIES[nz], marker=MARKERS[nz], ms=4, label=f"{nz} noise")
            ends.append([ys[-1], nz])
        ax.axhline(0.85, ls=(0, (4, 3)), lw=1, color=INK2)
        ax.text(1, 0.855, "G1b threshold", color=INK2, fontsize=8)
        ax.set_ylim(0.80, 1.02)
        # direct labels at line ends, spread apart so they never collide
        ends.sort(key=lambda e: e[0])
        gap = 0.012
        for k in range(1, len(ends)):
            ends[k][0] = max(ends[k][0], ends[k - 1][0] + gap)
        if ends and ends[-1][0] > 1.015:
            shift = ends[-1][0] - 1.015
            for e in ends:
                e[0] -= shift
        for y, nz in ends:
            ax.text(shares[-1] * 100 + 1.5, y, nz, color=INK, fontsize=8, va="center")
        ax.set_title(f"{reg} promises", fontsize=9, color=INK, loc="left")
        ax.set_xlabel("Capacity K/N (% of routes)")
        ax.set_xlim(0, 60)
        ax.grid(color=GRID, lw=0.6)
    axes[0][0].set_ylabel("Risk ranking / recoverability ranking")
    axes[0][-1].legend(frameon=False, fontsize=8, loc="center right", bbox_to_anchor=(1.0, 0.45))
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"gate1-capacity-curve{tag}.png"), dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("out_dir")
    ap.add_argument("--tag", default="")
    ap.add_argument("--figures", action="store_true")
    a = ap.parse_args()
    df = pd.read_csv(a.csv)
    os.makedirs(a.out_dir, exist_ok=True)
    summ = summarize(df)
    summ.to_csv(os.path.join(a.out_dir, f"summary{a.tag}.csv"), index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(summ.round(3).to_string(index=False))
    if {"base"} <= set(df.noise) and 0.35 in set(df.frac):
        outcome, lines = gate_decision(summ)
        print("\nGate 1 (pre-registered rule, f=0.35, base noise):")
        print("\n".join(lines))
        print("Outcome:", outcome)
    if a.figures:
        figures(df, a.out_dir, a.tag)


if __name__ == "__main__":
    main()
