"""Draw and freeze the Gate 3 route sample: whole station-days from the Amazon training split.

Station-days with at least 10 eligible routes (>= 10 stops) are shuffled with the seed and taken in
order until the sample holds at least 3,000 routes. Descriptive only: nothing is simulated.

Usage:
    python -I code/resequencing/gate3_sample.py <train_routes.pkl> <gate2c_main.csv> <out_csv>
"""
from __future__ import annotations

import pickle
import sys

import numpy as np
import pandas as pd

SEED, MIN_DAY, TARGET = 2029, 10, 3000


def main():
    src, g2c, out = sys.argv[1], sys.argv[2], sys.argv[3]
    routes = pickle.load(open(src, "rb"))
    df = pd.DataFrame([dict(ridx=i, route=r["id"], station=r["station"], date=r["date"], n=r["n"],
                            window_share=float(np.mean(~np.isnan(r["tw_close"][1:]))))
                       for i, r in enumerate(routes) if r["n"] - 1 >= 10])
    days = df.groupby(["station", "date"]).size()
    days = days[days >= MIN_DAY]
    order = np.random.default_rng(SEED).permutation(len(days))
    keep, total = [], 0
    for j in order:
        keep.append(days.index[j])
        total += int(days.iloc[j])
        if total >= TARGET:
            break
    keys = pd.MultiIndex.from_tuples(keep, names=["station", "date"])
    s = df.set_index(["station", "date"]).loc[keys].reset_index().sort_values("ridx")
    s.to_csv(out, index=False)
    prev = set(pd.read_csv(g2c, usecols=["route"]).route.unique())
    sizes = s.groupby(["station", "date"]).size()
    print(f"routes {len(s)}; station-days {len(sizes)}; stations {s.station.nunique()}; dates {s.date.nunique()}")
    print(f"routes per station-day: median {sizes.median():.0f}, min {sizes.min()}, max {sizes.max()}")
    print(f"stops per route: mean {s.n.mean() - 1:.1f}; stops with a window: {s.window_share.mean():.3f}")
    print(f"overlap with the Gate 2c sample: {s.route.isin(prev).sum()} routes")


if __name__ == "__main__":
    main()
