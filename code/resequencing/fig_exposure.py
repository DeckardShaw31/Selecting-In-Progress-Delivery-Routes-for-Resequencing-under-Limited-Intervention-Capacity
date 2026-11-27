"""Illustrative figure for the model section: two routes, equal remaining promise counts, Gaussian
deviations (sd 6*sqrt(i) min at the i-th remaining stop).

Route A: ETA promises, slack 60 min at every stop.
Route B: day-definite promises (slack 150 min) except two customer windows late in the sequence,
         at stops 10 and 11, with slack 0 and 10 min.

Usage:
    python -I code/resequencing/fig_exposure.py <out.png>   (also writes <out>.pdf, the vector copy used by the manuscript)
"""
import sys

import matplotlib
import numpy as np
from scipy.stats import norm

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"


def psi(z):
    return z * norm.cdf(z) + norm.pdf(z)


i = np.arange(1, 13)
sd = 6 * np.sqrt(i)
sA = np.full(12, 60.0)
sB = np.full(12, 150.0)
sB[9], sB[10] = 0.0, 10.0


def g(x, s):
    return np.sum(sd * psi((np.atleast_1d(x)[:, None] - s) / sd), axis=1)


def G(y, s):
    return np.sum(norm.cdf((np.atleast_1d(y)[:, None] - s) / sd), axis=1)


def main(out):
    plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.2))
    ax = axes[0]
    y = np.linspace(-40, 220, 600)
    ax.plot(y, G(y, sA), lw=2, color=BLUE, label="Route A: ETA promises")
    ax.plot(y, G(y, sB), lw=2, color=ORANGE, label="Route B: two late windows")
    ax.set_xlabel("Effective slack y = s − ξ (min)")
    ax.set_ylabel("Weighted distribution G_r(y) (promises)")
    ax.set_title("(a) Where the promises sit", fontsize=8.5, color=INK, loc="left")
    ax.grid(color=GRID, lw=0.6)
    ax = axes[1]
    x = np.linspace(-20, 70, 400)
    gA, gB = g(x, sA), g(x, sB)
    ax.plot(x, gA, lw=2, color=BLUE)
    ax.plot(x, gB, lw=2, color=ORANGE)
    DA, DB = 35.0, 0.0
    rA, rB = g(DA, sA)[0], g(DB, sB)[0]
    ax.plot([DA], [rA], "o", ms=6, color=BLUE, mec="white", mew=0.8, zorder=3)
    ax.plot([DB], [rB], "s", ms=6, color=ORANGE, mec="white", mew=0.8, zorder=3)
    ax.annotate(f"A: 35 min behind,\nrisk {rA:.1f} min", (DA, rA), xytext=(-6, 4), textcoords="offset points",
                fontsize=7.5, color=INK, ha="right", va="bottom")
    ax.annotate(f"B: on time,\nrisk {rB:.1f} min", (DB, rB), xytext=(-6, 10), textcoords="offset points",
                fontsize=7.5, color=INK, ha="right")
    ax.set_xlim(-20, 70)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Current delay x (min)")
    ax.set_ylabel("Exposure g_r(x) = expected lateness (min)")
    ax.set_title("(b) Exposure functions", fontsize=8.5, color=INK, loc="left")
    ax.grid(color=GRID, lw=0.6)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="upper center", ncol=2, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    plt.rcParams["pdf.fonttype"] = 42
    fig.savefig(out, dpi=600)
    fig.savefig(out.rsplit(".", 1)[0] + ".pdf")  # vector copy used by the manuscript
    print(f"risk A at 35 = {rA:.2f}; risk B at 0 = {rB:.2f}; crossing near x = {x[np.argmin(np.abs(gA - gB))]:.1f}")


if __name__ == "__main__":
    main(sys.argv[1])
