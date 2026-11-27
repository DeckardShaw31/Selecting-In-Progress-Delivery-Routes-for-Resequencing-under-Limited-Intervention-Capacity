"""Check 3c (tail form of Proposition 5) and bootstrap intervals for simulated efficiencies.

Implements analysis/theory/projection_tail_and_bootstrap.md (pre-specified 2026-11-27 before computing).
Reuses the route-level outputs of Gates 2c and 3; no re-simulation.

Usage:
    python -I code/resequencing/bootstrap_tail_check.py <project_root> <out_dir> [--B 2000]
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SHARE, SEED = 0.05, 2032


def topk_sum(target, score, K):
    return float(target[np.argsort(-score, kind="stable")[:K]].sum())


def K_of(N):
    return max(1, int(round(SHARE * N)))


def effs(c, idx=None):
    """Efficiencies of L1-L3 relative to L4 (equation 10) for cell arrays c, optionally on resampled rows idx."""
    get = (lambda k: c[k]) if idx is None else (lambda k: c[k][idx])
    d, rho = get("delta"), get("rho")
    K = K_of(len(d))
    v4 = topk_sum(d, rho, K)
    return {lv: topk_sum(d, get(lv), K) / v4 for lv in ("L1", "L2", "L3")}


# ------------------------------------------------------------------ Part 1: check 3c
def part1(root):
    rows = []
    for arm, f in (("A", "gate3_A.csv"), ("B", "gate3_B.csv")):
        df = pd.read_csv(os.path.join(root, "analysis/gate3/results", f))
        keys = ["sU", "regime", "frac", "noise"] if arm == "A" else ["noise", "phi", "regime"]
        for key, g in df.groupby(keys):
            g = g.sort_values("ridx")
            rho, pi, tau = g.rho_min.to_numpy(), g.L2_proj.to_numpy(), g.tau_ii_min.to_numpy()
            K = K_of(len(g))
            v_opt = topk_sum(rho, rho, K)
            rec = dict(zip(keys, key), arm=arm, N=len(g), K=K, V_opt=v_opt,
                       loss_L2=1 - topk_sum(rho, pi, K) / v_opt,
                       B_tail=float(np.sort(tau)[::-1][:2 * K].sum()) / v_opt,
                       B_max=2 * K * float(tau.max()) / v_opt,
                       omega=g.b_late_window_stops.sum() / g.b_late_stops.sum() if g.b_late_stops.sum() > 0 else np.nan)
            rows.append(rec)
    t = pd.DataFrame(rows)
    # cross-check against the stored Gate 3 cell tables
    for arm, f, keys in (("A", "gate3_cells_A.csv", ["sU", "regime", "frac", "noise"]),
                         ("B", "gate3_cells_B.csv", ["V_level", "phi", "regime"])):
        ref = pd.read_csv(os.path.join(root, "analysis/gate3", f))
        if arm == "B":
            ref = ref.rename(columns={"V_level": "noise"})
            keys = ["noise", "phi", "regime"]
        m = t[t.arm == arm].merge(ref[keys + ["loss_rho_L2", "B_ii"]], on=keys)
        assert len(m) == len(ref), (arm, len(m), len(ref))
        assert np.allclose(m.loss_L2, m.loss_rho_L2) and np.allclose(m.B_max, m.B_ii), f"mismatch with {f}"
    t["holds"] = t.B_tail >= t.loss_L2 - 1e-12
    t["below_max"] = t.B_tail <= t.B_max + 1e-12
    assert t.holds.all() and t.below_max.all(), "theorem violated: code error"
    sets = {"P (W cells, sU = 0)": t[(t.arm == "A") & (t.sU == 0)], "S1 (arm B)": t[t.arm == "B"],
            "S2 (arm A, sU > 0)": t[(t.arm == "A") & (t.sU > 0)]}
    summ = {}
    for name, s in sets.items():
        pos = s[s.loss_L2 > 0.01]
        summ[name] = dict(
            cells=len(s), spearman=float(spearmanr(s.B_tail, s.loss_L2).statistic),
            spearman_max_form=float(spearmanr(s.B_max, s.loss_L2).statistic),
            share_below_1=float((s.B_tail < 1).mean()),
            share_below_1_omega_ge_05=float((s[s.omega >= 0.5].B_tail < 1).mean()) if (s.omega >= 0.5).any() else None,
            share_below_1_omega_lt_05=float((s[s.omega < 0.5].B_tail < 1).mean()) if (s.omega < 0.5).any() else None,
            median_ratio_bound_to_loss=float(np.median(pos.B_tail / pos.loss_L2)) if len(pos) else None,
            B_tail_range=[float(s.B_tail.min()), float(s.B_tail.max())])
    return t, summ


# ------------------------------------------------------------------ Part 2: bootstrap
def cells_2c(root):
    df = pd.read_csv(os.path.join(root, "analysis/gate2c/results/gate2c_main.csv"),
                     usecols=["route", "ridx", "regime", "frac", "noise", "rho_min", "delta_min", "delay_plan_min",
                              "L2_proj", "L3_prior"])
    out, ids = {}, None
    for key, g in df.groupby(["regime", "frac", "noise"]):  # file order within the cell, as analyze_gate2c
        r = g.ridx.to_numpy()
        ids = np.sort(r) if ids is None else ids
        assert np.array_equal(np.sort(r), ids), "route sets differ across cells"
        pos = np.empty(ids.max() + 1, int)
        pos[r] = np.arange(len(r))  # ridx -> row position in this cell's file order
        out[key] = dict(delta=g.delta_min.to_numpy(), rho=g.rho_min.to_numpy(), L1=g.delay_plan_min.to_numpy(),
                        L2=g.L2_proj.to_numpy(), L3=g.L3_prior.to_numpy(), pos=pos, ids=ids)
    return out


def cells_3(root):
    a = pd.read_csv(os.path.join(root, "analysis/gate3/results/gate3_A.csv"))
    b = pd.read_csv(os.path.join(root, "analysis/gate3/results/gate3_B.csv"))
    out = {}
    sel = [(("A", 0.5, r, 0.35, nz), a[(a.sU == 0.5) & (a.regime == r) & (a.frac == 0.35) & (a.noise == nz)])
           for r in a.regime.unique() for nz in ("base", "high")]
    sel.append((("B", "base", 0.75, "MIXB"), b[(b.noise == "base") & np.isclose(b.phi, 0.75) & (b.regime == "MIXB")]))
    for key, g in sel:
        if len(g) == 0:
            continue
        g = g.sort_values("ridx")  # as analyze_gate3.cell_record
        out[key] = dict(delta=g.delta_min.to_numpy(), rho=g.rho_min.to_numpy(), L1=g.delay_plan_min.to_numpy(),
                        L2=g.L2_proj.to_numpy(), L3=g.L3_prior.to_numpy(),
                        cluster=(g.station.astype(str) + "|" + g.date.astype(str)).to_numpy())
    return out


def ci(x):
    return [float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))]


def bootstrap(cells, B, cluster):
    rng = np.random.default_rng(SEED)
    first = next(iter(cells.values()))
    N = len(first["delta"])
    if cluster:
        labs = first["cluster"]
        for c in cells.values():
            assert np.array_equal(c["cluster"], labs), "cluster labels differ across cells"
        uniq, inv = np.unique(labs, return_inverse=True)
        members = [np.flatnonzero(inv == k) for k in range(len(uniq))]
    draws = {k: {lv: [] for lv in ("L1", "L2", "L3", "L2-L1")} for k in cells}
    for _ in range(B):
        if cluster:
            pick = rng.integers(0, len(members), len(members))
            idx = np.sort(np.concatenate([members[k] for k in pick]))
        else:  # one resample of routes for all cells; rows kept in each cell's own order (ties as in the estimate)
            draw = first["ids"][rng.integers(0, N, N)]
        for k, c in cells.items():
            e = effs(c, idx if cluster else np.sort(c["pos"][draw]))
            for lv in ("L1", "L2", "L3"):
                draws[k][lv].append(e[lv])
            draws[k]["L2-L1"].append(e["L2"] - e["L1"])
    rows = []
    for k, c in cells.items():
        e0 = effs(c)
        rec = dict(cell=" / ".join(map(str, k)))
        for lv in ("L1", "L2", "L3"):
            rec[f"{lv}"] = e0[lv]
            rec[f"{lv}_lo"], rec[f"{lv}_hi"] = ci(draws[k][lv])
        rec["L2-L1"] = e0["L2"] - e0["L1"]
        rec["L2-L1_lo"], rec["L2-L1_hi"] = ci(draws[k]["L2-L1"])
        rows.append(rec)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("out_dir")
    ap.add_argument("--B", type=int, default=2000)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    pd.set_option("display.width", 250)

    t, summ = part1(a.root)
    t.to_csv(os.path.join(a.out_dir, "check3c_cells.csv"), index=False)
    print("Part 1 (check 3c): all bounds hold:", bool(t.holds.all()), "| tail <= max form:", bool(t.below_max.all()))
    print(json.dumps(summ, indent=1))

    c2 = cells_2c(a.root)
    ref = pd.read_csv(os.path.join(a.root, "analysis/gate2c/gate2c_cells.csv"))
    for (reg, f, nz), c in c2.items():
        r = ref[(ref.regime == reg) & np.isclose(ref.frac, f) & (ref.noise == nz)].iloc[0]
        e = effs(c)
        assert all(abs(e[lv] - r[f"{lv}_eff_5"]) < 1e-12 for lv in ("L1", "L2", "L3")), ("2c point estimate", reg, f, nz)
    b2 = bootstrap(c2, a.B, cluster=False)
    b2.to_csv(os.path.join(a.out_dir, "bootstrap_gate2c.csv"), index=False)
    c3 = cells_3(a.root)
    refA = pd.read_csv(os.path.join(a.root, "analysis/gate3/gate3_cells_A.csv"))
    for k, c in c3.items():
        if k[0] == "A":
            r = refA[np.isclose(refA.sU, 0.5) & (refA.regime == k[2]) & np.isclose(refA.frac, 0.35) & (refA.noise == k[4])].iloc[0]
            e = effs(c)
            assert all(abs(e[lv] - r[f"{lv}_eff_5"]) < 1e-12 for lv in ("L1", "L2", "L3")), ("3 point estimate", k)
    b3 = bootstrap(c3, a.B, cluster=True)
    b3.to_csv(os.path.join(a.out_dir, "bootstrap_gate3.csv"), index=False)
    print("\nPart 2, Gate 2c (route bootstrap, B = %d):" % a.B)
    print(b2.round(3).to_string(index=False))
    print("\nPart 2, Gate 3 (station-day cluster bootstrap):")
    print(b3.round(3).to_string(index=False))
    hw = {lv: float(((b2[f"{lv}_hi"] - b2[f"{lv}_lo"]) / 2).max()) for lv in ("L1", "L2", "L3", "L2-L1")}
    json.dump(dict(check3c=summ, gate2c_max_halfwidth=hw, B=a.B, seed=SEED),
              open(os.path.join(a.out_dir, "summary.json"), "w"), indent=1)
    print("\nlargest 95% half-width across Gate 2c cells:", {k: round(v, 3) for k, v in hw.items()})


if __name__ == "__main__":
    main()
