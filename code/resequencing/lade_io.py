"""Load and clean LaDe pickup tasks for Gate 4 (see analysis/gate4/preregistration.md).

Times are parsed as minutes from 2020-01-01 (a leap year, so 02-29 parses; the data span about
six months within one year). A courier-day (city, courier_id, ds) is excluded whole if any of its
tasks has a missing or unparseable time or coordinate, window end < window start, pickup < accept,
or lateness > 24 h. Eligible courier-days have >= 10 tasks after cleaning.

Usage:
    python -I code/resequencing/lade_io.py <raw_pickup_dir> <out.pkl> <out_cleaning.json>
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

TIME_COLS = ["accept_time", "time_window_start", "time_window_end", "pickup_time"]
REF = pd.Timestamp("2020-01-01")


def to_min(s):
    t = pd.to_datetime("2020-" + s.astype(str), format="%Y-%m-%d %H:%M:%S", errors="coerce")
    return (t - REF).dt.total_seconds() / 60.0


def main():
    raw, out, out_json = sys.argv[1], sys.argv[2], sys.argv[3]
    parts = []
    for f in sorted(glob.glob(os.path.join(raw, "pickup_*.csv"))):
        d = pd.read_csv(f, usecols=["order_id", "city", "courier_id", "lng", "lat", "ds"] + TIME_COLS)
        parts.append(d)
    d = pd.concat(parts, ignore_index=True)
    rep = dict(tasks_raw=int(len(d)))
    for c in TIME_COLS:
        d[c] = to_min(d[c])
    d["row"] = np.arange(len(d))
    key = ["city", "courier_id", "ds"]
    rep["courier_days_raw"] = int(d.groupby(key).ngroups)
    bad = {
        "missing_or_unparseable": d[TIME_COLS + ["lng", "lat"]].isna().any(axis=1) | d.courier_id.isna() | d.ds.isna(),
        "window_end_before_start": d.time_window_end < d.time_window_start,
        "pickup_before_accept": d.pickup_time < d.accept_time,
        "late_over_24h": (d.pickup_time - d.time_window_end) > 24 * 60,
    }
    excl = pd.Series(False, index=d.index)
    for name, m in bad.items():
        m = m.fillna(False)
        days = d.loc[m, key].drop_duplicates()
        rep[f"courier_days_with_{name}"] = int(len(days))
        rep[f"tasks_flagged_{name}"] = int(m.sum())
        excl |= m
    bad_days = d.loc[excl, key].drop_duplicates()
    d = d.merge(bad_days.assign(_bad=1), on=key, how="left")
    d = d[d._bad.isna()].drop(columns="_bad")
    rep["courier_days_after_cleaning"] = int(d.groupby(key).ngroups)
    size = d.groupby(key)["order_id"].transform("size")
    d = d[size >= 10].copy()
    rep["courier_days_eligible"] = int(d.groupby(key).ngroups)
    rep["tasks_eligible"] = int(len(d))
    # split by date within city: first half of distinct ds = calibration
    d["split"] = "evaluation"
    for city, g in d.groupby("city"):
        ds = sorted(g.ds.unique())
        cal = set(ds[: len(ds) // 2])
        d.loc[(d.city == city) & d.ds.isin(cal), "split"] = "calibration"
        rep[f"dates_{city}"] = [int(ds[0]), int(ds[-1]), len(ds)]
    rep["courier_days_by_split"] = d.groupby("split").apply(lambda g: g.groupby(key).ngroups).astype(int).to_dict()
    d = d.sort_values(key + ["pickup_time", "row"]).reset_index(drop=True)
    d.to_pickle(out)
    json.dump(rep, open(out_json, "w"), indent=1)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
