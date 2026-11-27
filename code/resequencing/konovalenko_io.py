"""Build the Konovalenko et al. (2026) stop table used by Gate 2d and the Gate 3 calibration.

Source: Konovalenko, Hvattum & Iversen (2026), Mendeley Data, doi 10.17632/kkwgfvmtxn.1, file
routes_performance.xlsx (CC BY 4.0). The table is stored as read, one row per stop visit, with no cleaning:
eligibility rules (>= 10 non-depot stops, complete times), the calibration/evaluation split by week, and the
pace fit are applied downstream in gate2d.py and gate3_calibrate.py, as pre-specified.

Columns: Route ID, Driver ID, Stop ID, Address ID, Week ID, Country, Day of Week, IndexP (planned position),
IndexA (actual position), Arrived Time, Earliest Time, Latest Time (minutes), DistanceP, DistanceA (km to the
previous stop in the planned and actual order), Depot, Delivery.

Usage:
    python -I code/resequencing/konovalenko_io.py <routes_performance.xlsx> <out.pkl>

Check (2026-11-27): the output equals data/processed/konovalenko_stops.pkl used for the paper
(pandas.testing.assert_frame_equal).
"""
from __future__ import annotations

import sys

import pandas as pd

COLUMNS = ["Route ID", "Driver ID", "Stop ID", "Address ID", "Week ID", "Country", "Day of Week", "IndexP", "IndexA",
           "Arrived Time", "Earliest Time", "Latest Time", "DistanceP", "DistanceA", "Depot", "Delivery"]


def main():
    src, out = sys.argv[1], sys.argv[2]
    d = pd.read_excel(src, sheet_name="Sheet1")
    assert list(d.columns) == COLUMNS, f"unexpected columns: {list(d.columns)}"
    d.to_pickle(out)
    nd = d[d["Depot"] == 0]
    late = (nd["Arrived Time"] > nd["Latest Time"]).mean()
    print(f"{len(d):,} stop visits, {d['Route ID'].nunique():,} routes, {d['Driver ID'].nunique()} drivers, "
          f"{d['Week ID'].nunique()} weeks; non-depot visits after the latest allowed time: {late:.3f}")
    if len(sys.argv) > 3:  # optional: compare with a reference pickle
        pd.testing.assert_frame_equal(d, pd.read_pickle(sys.argv[3]))
        print("identical to", sys.argv[3])


if __name__ == "__main__":
    main()
