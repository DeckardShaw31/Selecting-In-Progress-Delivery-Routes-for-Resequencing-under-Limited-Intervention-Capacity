"""Gate 3: imperfect prediction, delay structure and the loss bound (see analysis/gate3/preregistration.md).

Arm A (E6): the Amazon means are the carrier's predictions (plan, promises, L1 and all forecasts);
true means carry persistent multiplicative errors with log-SD sigma_U. The oracle, the risk L4*
and the value use the true law. At sigma_U = 0 (W) the dispatcher's dispersion parameters are also
scaled by 0.5 and 2 (misspecified dispersion).
Arm B (E7): perfectly specified; travel noise split into a route-wide share phi of a fixed total
log-variance V (base or high level).

Usage:
    python -I code/resequencing/gate3.py <train_routes.pkl> <sample_routes.csv> <out.csv> --arm A|B
           [--workers W] [--seed 2029] [--limit N]
"""
from __future__ import annotations

import argparse
import os
import pickle
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate1 import (ETA_HALF_WIDTH, MAX_PASSES, NOISE, S_IN, S_OOS, nominal_eta,  # noqa: E402
                   nominal_schedule_wait, promise_arrays, regime_windows)
from sim import draw_future, eval_seq, local_search, posterior_theta, realize_route  # noqa: E402

NOISE_IDX = {k: i for i, k in enumerate(NOISE)}
REG_W = ["ETA", "MIXA", "DAY", "MIXB", "MIXB20", "MIXB50", "TW"]
REG_U = ["ETA", "DAY", "MIXB", "MIXB20", "MIXB50", "TW"]
SIGMA_U = [0.25, 0.50, 1.00]
REG_B = ["DAY", "MIXB", "MIXB50", "TW"]
PHI = [0.10, 0.25, 0.50, 0.75, 0.90]
MISSPEC = [0.5, 2.0]

COLS = ["route", "ridx", "station", "date", "n", "U", "arm", "noise", "sU", "phi", "sR", "sE", "sS", "frac", "regime",
        "t0_min", "delay_so_far_min", "delay_plan_min", "theta_true", "theta_post_true", "theta_post_carrier",
        "rho_min", "lam_min", "delta_min", "rho_in_min", "best_in_min", "passes", "stops_moved",
        "L2_proj", "L3_prior", "L4f", "b_late_stops", "b_late_window_stops", "windows_remaining",
        "window_share_route", "tau_ii_min", "L3_m05", "L4_m05", "L3_m2", "L4_m2", "secs"]

_ROUTES = None


def _init(path):
    global _ROUTES
    with open(path, "rb") as fh:
        _ROUTES = pickle.load(fh)


@njit(cache=True)
def start_times(seq, t0, T, Sv, opn, wait):
    """Service-start times (after any waiting) of the stops in `seq`, per scenario: (S, len(seq))."""
    S = T.shape[0]
    U = seq.shape[0]
    out = np.empty((S, U))
    for s in range(S):
        t = t0
        prev = 0
        for k in range(U):
            j = seq[k]
            t += T[s, prev, j]
            if wait and t < opn[j]:
                t = opn[j]
            out[s, k] = t
            t += Sv[s, j]
            prev = j
    return out


def true_means(r, sU, seed, ridx):
    """True mean travel and service times: carrier means x exp(sU z - sU^2/2), with z drawn once per
    route (the same z for every sU, regime, check-in and noise level)."""
    tb = r["tbar"].astype(np.float64)
    sb = r["sbar"].astype(np.float64)
    if sU == 0:
        return tb, sb
    rng = np.random.default_rng([seed, ridx, 7])
    n = r["n"]
    Z = rng.normal(size=(n, n))
    Zs = rng.normal(size=n)
    tb_t = tb * np.exp(sU * Z - sU ** 2 / 2)
    sb_t = sb * np.exp(sU * Zs - sU ** 2 / 2)
    sb_t[0] = 0.0
    return tb_t, sb_t


def _obs(legs, tb, m, corr):
    return np.array([np.log(legs[k] / tb[k - 1, k]) + corr for k in range(1, m + 1) if tb[k - 1, k] > 0])


def run_config(t):
    ridx, seed, lv = t["ridx"], t["seed"], t["rng_level"]
    sR, sE, sS, sU = t["sR"], t["sE"], t["sS"], t["sU"]
    r = _ROUTES[ridx]
    n = r["n"]
    eta_hw = ETA_HALF_WIDTH
    tb_c = r["tbar"].astype(np.float64)
    sb_c = r["sbar"]
    tb_t, sb_t = true_means(r, sU, seed, ridx)
    legs, svc, theta_true, _ = realize_route(np.random.default_rng([seed, ridx, 0, lv]), tb_t, sb_t, sR, sE, sS)
    eta = nominal_eta(r)  # carrier plan
    WIN = {}
    for regime in t["regimes"]:
        o_full, c_full = regime_windows(r, regime, eta, seed, ridx)
        st_w, en_w = nominal_schedule_wait(r, o_full)
        WIN[regime] = (o_full, c_full, st_w, en_w + eta_hw)
    corr = (sR ** 2 + sE ** 2) / 2.0
    out = []
    for f in t["fracs"]:
        m = int(round(f * (n - 1)))
        if m < 1 or n - 1 - m < 2:
            continue
        idx = np.arange(m, n)
        U = len(idx) - 1
        loc = np.ix_(idx, idx)
        tl_t, sl_t = tb_t[loc], sb_t[idx].copy()
        tl_c, sl_c = tb_c[loc], sb_c[idx].copy()
        mu_t, var_t = posterior_theta(_obs(legs, tb_t, m, corr), sR, sE)
        mu_c, var_c = posterior_theta(_obs(legs, tb_c, m, corr), sR, sE)
        fk = int(f * 100)

        def rng(stream):
            return np.random.default_rng([seed, ridx, stream, lv, fk])

        # regime-specific promises and realized time at check-in (true execution, carrier promises)
        PR = {}
        cur = np.arange(1, U + 1, dtype=np.int64)
        for regime in t["regimes"]:
            o_full, c_full, start_w, day_deadline = WIN[regime]
            opn, cls, w, wait = promise_arrays(r, idx, regime, eta, eta_hw, day_deadline, start_w, o_full, c_full)
            t0 = 0.0
            for k in range(1, m + 1):
                t0 += legs[k]
                if wait and not np.isnan(o_full[k]) and t0 < o_full[k]:
                    t0 = o_full[k]
                t0 += svc[k]
            PR[regime] = (opn, cls, w, wait, t0, o_full, c_full, start_w)
        rows = {}
        # --- truth: oracle resequencing and risk (in-sample and out-of-sample under the true law)
        T_in, Sv_in = draw_future(rng(1), tl_t, sl_t, mu_t, var_t, sR, sE, sS, S_IN)
        T_oos, Sv_oos = draw_future(rng(2), tl_t, sl_t, mu_t, var_t, sR, sE, sS, S_OOS)
        T_det, Sv_det = tl_c[None, :, :].copy(), sl_c[None, :].copy()
        Sv_det[:, 0] = 0.0
        for regime in t["regimes"]:
            opn, cls, w, wait, t0, o_full, c_full, start_w = PR[regime]
            tic = time.perf_counter()
            rho_in = eval_seq(cur, t0, T_in, Sv_in, opn, cls, w, wait, 0.0)
            best, best_in, passes = cur, rho_in, 0
            if rho_in > 0:
                s1, c1, p1 = local_search(cur, t0, T_in, Sv_in, opn, cls, w, wait, MAX_PASSES, 0.0)
                best, best_in, passes = s1, c1, p1
                if regime != "ETA":
                    edd = (np.argsort(cls[1:], kind="stable") + 1).astype(np.int64)
                    if not np.array_equal(edd, cur):
                        s2, c2, p2 = local_search(edd, t0, T_in, Sv_in, opn, cls, w, wait, MAX_PASSES, 0.0)
                        if c2 < best_in:
                            best, best_in, passes = s2, c2, passes + p2
            rho = eval_seq(cur, t0, T_oos, Sv_oos, opn, cls, w, wait, 0.0)
            lam = eval_seq(best, t0, T_oos, Sv_oos, opn, cls, w, wait, 0.0) if passes > 0 else rho
            b_cnt = eval_seq(cur, t0, T_oos, Sv_oos, opn, cls, w, wait, 1.0)
            has_win = ~np.isnan(c_full[idx])
            has_win[0] = False
            if regime == "ETA":
                has_win[:] = False
            b_win = eval_seq(cur, t0, T_oos, Sv_oos, opn, cls, w * has_win, wait, 1.0) if has_win.any() else 0.0
            L2 = eval_seq(cur, t0, T_det, Sv_det, opn, cls, w, wait, 0.0)
            # Proposition 2(ii) term: sum_i w_i E|A_i - Ahat_i| 1{max(A_i, Ahat_i) > d_i}
            A = start_times(cur, t0, T_oos, Sv_oos, opn, wait)
            Ah = start_times(cur, t0, T_det, Sv_det, opn, wait)[0]
            d = cls[cur]
            ind = np.maximum(A, Ah[None, :]) > d[None, :]
            tau = float(np.sum(w[cur] * np.mean(np.abs(A - Ah[None, :]) * ind, axis=0)))
            plan_done = (start_w[m] if wait else eta[m]) + sb_c[m]
            rows[regime] = dict(
                route=r["id"], ridx=ridx, station=r["station"], date=r["date"], n=n, U=U, arm=t["arm"],
                noise=t["noise"], sU=sU, phi=t["phi"], sR=sR, sE=sE, sS=sS, frac=f, regime=regime,
                t0_min=t0 / 60, delay_so_far_min=(t0 - (eta[m] + sb_c[m])) / 60, delay_plan_min=(t0 - plan_done) / 60,
                theta_true=theta_true, theta_post_true=mu_t, theta_post_carrier=mu_c,
                rho_min=rho / 60, lam_min=lam / 60, delta_min=(rho - lam) / 60, rho_in_min=rho_in / 60,
                best_in_min=best_in / 60, passes=passes, stops_moved=int(np.sum(best != cur)),
                L2_proj=L2 / 60, L4f=rho / 60, b_late_stops=b_cnt, b_late_window_stops=b_win,
                windows_remaining=int(np.sum(has_win)), window_share_route=float(np.mean(~np.isnan(c_full[1:]))),
                tau_ii_min=tau / 60, L3_m05=np.nan, L4_m05=np.nan, L3_m2=np.nan, L4_m2=np.nan,
                secs=time.perf_counter() - tic,
            )
        del T_in, Sv_in, T_oos, Sv_oos
        # --- dispatcher: prior forecast with carrier means
        T_x, Sv_x = draw_future(rng(3), tl_c, sl_c, 0.0, sR ** 2, sR, sE, sS, S_OOS)
        for regime in t["regimes"]:
            opn, cls, w, wait, t0 = PR[regime][:5]
            rows[regime]["L3_prior"] = eval_seq(cur, t0, T_x, Sv_x, opn, cls, w, wait, 0.0) / 60
        # --- dispatcher: posterior forecast with carrier means (equals L4* when sU = 0)
        if sU > 0:
            T_x, Sv_x = draw_future(rng(2), tl_c, sl_c, mu_c, var_c, sR, sE, sS, S_OOS)
            for regime in t["regimes"]:
                opn, cls, w, wait, t0 = PR[regime][:5]
                rows[regime]["L4f"] = eval_seq(cur, t0, T_x, Sv_x, opn, cls, w, wait, 0.0) / 60
        # --- misspecified dispersion (W only): all sigmas x c in the dispatcher's model
        if t["misspec"]:
            for i, c in enumerate(MISSPEC):
                tag = "m05" if c == 0.5 else "m2"
                cR, cE, cS = c * sR, c * sE, c * sS
                mu_m, var_m = posterior_theta(_obs(legs, tb_c, m, (cR ** 2 + cE ** 2) / 2.0), cR, cE)
                T_x, Sv_x = draw_future(rng(4 + i), tl_c, sl_c, 0.0, cR ** 2, cR, cE, cS, S_OOS)
                for regime in t["regimes"]:
                    opn, cls, w, wait, t0 = PR[regime][:5]
                    rows[regime][f"L3_{tag}"] = eval_seq(cur, t0, T_x, Sv_x, opn, cls, w, wait, 0.0) / 60
                T_x, Sv_x = draw_future(rng(6 + i), tl_c, sl_c, mu_m, var_m, cR, cE, cS, S_OOS)
                for regime in t["regimes"]:
                    opn, cls, w, wait, t0 = PR[regime][:5]
                    rows[regime][f"L4_{tag}"] = eval_seq(cur, t0, T_x, Sv_x, opn, cls, w, wait, 0.0) / 60
        del T_x, Sv_x
        out += [rows[g] for g in t["regimes"]]
    return out


def tasks_for(ridx_list, arm, seed):
    tasks = []
    for ridx in ridx_list:
        if arm == "A":
            for nz in ("base", "high"):
                sR, sE, sS = NOISE[nz]
                base = dict(ridx=ridx, seed=seed, arm="A", noise=nz, rng_level=NOISE_IDX[nz], sR=sR, sE=sE, sS=sS,
                            phi=sR ** 2 / (sR ** 2 + sE ** 2))
                tasks.append(dict(base, sU=0.0, fracs=[0.35, 0.50], regimes=REG_W, misspec=True))
                for sU in SIGMA_U:
                    tasks.append(dict(base, sU=sU, fracs=[0.35], regimes=REG_U, misspec=False))
        else:
            for j, lvl in enumerate(("base", "high")):
                sR0, sE0, sS0 = NOISE[lvl]
                V = sR0 ** 2 + sE0 ** 2
                for phi in PHI:
                    tasks.append(dict(ridx=ridx, seed=seed, arm="B", noise=lvl, rng_level=20 + j, sU=0.0, phi=phi,
                                      sR=float(np.sqrt(phi * V)), sE=float(np.sqrt((1 - phi) * V)), sS=sS0,
                                      fracs=[0.35], regimes=REG_B, misspec=False))
    return tasks


def _key(d):
    return (int(d["ridx"]), d["noise"], round(float(d["sU"]), 3), round(float(d["phi"]), 3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("routes")
    ap.add_argument("sample")
    ap.add_argument("out")
    ap.add_argument("--arm", required=True, choices=["A", "B"])
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--seed", type=int, default=2029)
    ap.add_argument("--limit", type=int, default=0, help="first N sample routes only (testing)")
    a = ap.parse_args()
    ridx = pd.read_csv(a.sample).ridx.astype(int).tolist()
    if a.limit:
        ridx = ridx[: a.limit]
    tasks = tasks_for(ridx, a.arm, a.seed)
    done = set()
    if os.path.exists(a.out):
        prev = pd.read_csv(a.out, usecols=["ridx", "noise", "sU", "phi"], on_bad_lines="skip")
        done = {_key(row) for row in prev.drop_duplicates().to_dict("records")}
        tasks = [t for t in tasks if _key(t) not in done]
    else:
        pd.DataFrame(columns=COLS).to_csv(a.out, index=False)
    print(f"arm {a.arm}; routes {len(ridx)}; tasks to run {len(tasks)} (done {len(done)}); workers {a.workers}", flush=True)
    t0 = time.time()
    buf = []
    with Pool(a.workers, initializer=_init, initargs=(a.routes,)) as pool:
        for k, res in enumerate(pool.imap_unordered(run_config, tasks, chunksize=1)):
            buf += res
            if (k + 1) % 200 == 0 or k + 1 == len(tasks):
                pd.DataFrame(buf, columns=COLS).to_csv(a.out, mode="a", header=False, index=False)
                buf = []
                el = time.time() - t0
                print(f"{k + 1}/{len(tasks)} tasks, {el:.0f}s, eta {el / (k + 1) * (len(tasks) - k - 1) / 60:.0f} min", flush=True)
    print(f"done in {time.time() - t0:.0f}s -> {a.out}")


if __name__ == "__main__":
    main()
