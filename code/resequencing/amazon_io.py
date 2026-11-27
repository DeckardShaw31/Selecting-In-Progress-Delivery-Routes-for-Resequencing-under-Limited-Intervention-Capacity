"""Load the Amazon Last Mile Routing Research Challenge data into compact per-route arrays.

Usage:
    python -I code/resequencing/amazon_io.py <raw_dir> <out_pickle> [--prefix eval_]

Each route is stored in *actual sequence order* (position 0 = station):
    id, station, date, n, stop_ids, lat, lng, zone,
    tbar  (n, n) float32  mean travel time between positions (seconds)        [O]
    sbar  (n,)   float64  stop service time = sum of package planned times    [D]
    tw_open, tw_close (n,) seconds after departure; NaN if no window          [D]
    n_pkg (n,)   packages per stop
Time windows: stop open = latest package window start, close = earliest window end.
"""
from __future__ import annotations

import json
import mmap
import math
import pickle
import re
import sys
from datetime import datetime, timezone

import numpy as np


def _ts(s):
    if s is None or (isinstance(s, float) and math.isnan(s)):
        return None
    return datetime.strptime(s, "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc).timestamp()


def iter_top_level(path, key_prefix=b'"RouteID_'):
    """Yield (route_id, value_dict) from a large {RouteID: {...}} JSON file without loading it whole.

    Relies on RouteID_ appearing only as top-level keys (true for the travel-time file,
    whose inner keys are two-letter stop codes).
    """
    with open(path, "rb") as fh:
        mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
        offs = [m.start() for m in re.finditer(re.escape(key_prefix), mm)]
        for k, o in enumerate(offs):
            end = offs[k + 1] if k + 1 < len(offs) else len(mm)
            chunk = mm[o:end].decode("utf-8").rstrip()
            if k + 1 == len(offs):
                chunk = chunk[:-1].rstrip()  # drop top-level closing brace
            chunk = chunk.rstrip(",").rstrip()
            d = json.loads("{" + chunk + "}")
            (rid, val), = d.items()
            yield rid, val
        mm.close()


def build(raw_dir, out_path, prefix="eval_"):
    route_data = json.load(open(f"{raw_dir}/{prefix}route_data.json", encoding="utf-8"))
    package_data = json.load(open(f"{raw_dir}/{prefix}package_data.json", encoding="utf-8"))
    actual = json.load(open(f"{raw_dir}/{prefix}actual_sequences.json", encoding="utf-8"))
    routes = []
    skipped = {"no_actual": 0, "bad_seq": 0}
    for rid, tt in iter_top_level(f"{raw_dir}/{prefix}travel_times.json"):
        rd = route_data[rid]
        if rid not in actual:
            skipped["no_actual"] += 1
            continue
        order = actual[rid]["actual"]
        stop_ids = sorted(order, key=lambda s: order[s])
        n = len(stop_ids)
        if sorted(order.values()) != list(range(n)) or rd["stops"][stop_ids[0]]["type"] != "Station":
            skipped["bad_seq"] += 1
            continue
        dep = datetime.strptime(f'{rd["date_YYYY_MM_DD"]} {rd["departure_time_utc"]}', "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc).timestamp()
        tbar = np.zeros((n, n), dtype=np.float32)
        for a, sa in enumerate(stop_ids):
            row = tt[sa]
            for b, sb in enumerate(stop_ids):
                tbar[a, b] = row[sb]
        sbar = np.zeros(n)
        npk = np.zeros(n, dtype=np.int32)
        opn = np.full(n, np.nan)
        cls = np.full(n, np.nan)
        pk = package_data.get(rid, {})
        for a, s in enumerate(stop_ids):
            for p in pk.get(s, {}).values():
                sbar[a] += float(p.get("planned_service_time_seconds") or 0.0)
                npk[a] += 1
                tw = p.get("time_window") or {}
                st, en = _ts(tw.get("start_time_utc")), _ts(tw.get("end_time_utc"))
                if st is not None:
                    opn[a] = st - dep if np.isnan(opn[a]) else max(opn[a], st - dep)
                if en is not None:
                    cls[a] = en - dep if np.isnan(cls[a]) else min(cls[a], en - dep)
        routes.append(dict(
            id=rid, station=rd["station_code"], date=rd["date_YYYY_MM_DD"], dep_utc=dep, n=n,
            stop_ids=stop_ids,
            lat=np.array([rd["stops"][s]["lat"] for s in stop_ids]),
            lng=np.array([rd["stops"][s]["lng"] for s in stop_ids]),
            zone=[rd["stops"][s].get("zone_id") for s in stop_ids],
            tbar=tbar, sbar=sbar, tw_open=opn, tw_close=cls, n_pkg=npk,
        ))
    with open(out_path, "wb") as fh:
        pickle.dump(routes, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"routes: {len(routes)}; skipped: {skipped}")


if __name__ == "__main__":
    raw, out = sys.argv[1], sys.argv[2]
    prefix = sys.argv[4] if len(sys.argv) > 4 and sys.argv[3] == "--prefix" else "eval_"
    build(raw, out, prefix)
