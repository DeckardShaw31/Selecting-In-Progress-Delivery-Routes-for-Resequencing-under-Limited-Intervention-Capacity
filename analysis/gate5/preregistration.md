# Gate 5 pre-registration: delay × delay sensitivity

Written 2026-11-26, before any Gate 5 score was computed on any route.

## Purpose

Lemma 3 (manuscript Section 3.6) shows that at large delays risk behaves like $\tilde W_r(D_r-\tilde\mu_r)$. Here $\tilde W_r=\sum_i w_{ri}\lambda_{ri}$ is the route's **delay sensitivity**: the risk added per minute of delay once its remaining promises are broken. Behind-schedule triage ranks by $D_r$ alone and ignores $\tilde W_r$. By Proposition 3, with unequal sensitivities no finite guarantee holds for delay ranking.

The suggested rule ranks routes by **delay × delay sensitivity**. It needs no promise positions or slacks, and in its simplest form no forecasting model. Gate 5 asks:
- how much of behind-schedule triage's shortfall the rule recovers;
- whether the shortfall comes from unequal sensitivities (which the rule fixes) or from heterogeneous slack (which it cannot fix, and which the projection fixes).

## What has and has not been seen

- **Already known:** the efficiencies of L0–L4 in every simulated cell, and the observed-data results of B1, B2 and F.
- **Never computed:** the three Gate 5 scores, on any route.
- **Data.** The test uses only existing simulation outputs and route data. Nothing is re-simulated, so the value of every selection (out-of-sample removable lateness `delta_min`) is fixed.

## Scores (computed per route and cell from information available to a dispatcher)

- **D:** delay against the plan (`delay_plan_min`); this is L1.
- **W_r:** remaining promise weight, the sum of the promise weights over the remaining stops after the current one.
  - Every remaining stop under ETA, MIXA, DAY, MIXB, MIXB20 and MIXB50 (`U`).
  - Remaining stops with a customer window under TW (`windows_remaining`).
- **S1, delay × remaining promises (primary):** W_r · D_r. It is model-free and uses only the count of remaining promises.
- **S2, pace-extrapolated delay (secondary):** (D_r / Λ_r) · Σ_i w_ri (Λ_r + c_ri). This is the perfect-learning sensitivity of Proposition 7(c).
  - Λ_r is the planned travel time of the completed legs.
  - c_ri is the planned travel time from the check-in to remaining stop i, using the carrier's mean travel times from the route data.
  - If Λ_r = 0, S2 = S1.
- **S3, model-based sensitivity (secondary):** D_r · Σ_i w_ri (1 + c_ri γ_r), with γ_r from Proposition 7.
  - γ_r = Λ_r φ σ² / Var D_r and Var D_r = σ²[Λ_r² φ + Q_r(1 − φ)] + σ_S² Σ_{k ≤ m} s̄_k², in seconds.
  - σ² = σ_R² + σ_E² and φ = σ_R²/σ² are the cell's (dispatcher's) noise parameters.
  - Q_r is the sum of squared completed planned leg times, and s̄_k are the planned service times of the completed stops.

**Value.** As in Gates 2c and 3: the removable lateness (`delta_min`) of the K routes ranked highest by the score, divided by that of the K routes ranked highest by risk (L4). The tie rule is a stable sort by route index. The primary capacity is K/N = 5%; 2% and 10% are reported.

**Gap recovered.** $G(S)=[\mathrm{eff}(S)-\mathrm{eff}(L1)]/[1-\mathrm{eff}(L1)]$, computed in cells where $1-\mathrm{eff}(L1)\ge0.02$.

## Data sets

| Set | Cells | Role |
|---|---|---|
| Gate 2c confirmatory sample (2,999 routes; 7 regimes × 2 check-ins × 2 noise levels) | 28 | **Primary** |
| Gate 3 W (3,028 routes, same regimes) | 28 | Replication |
| Gate 3 Arm B (delay structure, φ grid) | 40 | Secondary, no threshold |
| Gate 3 σ_U > 0 (imperfect prediction) | 36 | Secondary, no threshold |
| Gate 4 LaDe pickups (evaluation courier-days, check-in 35%) | 1 | Observed, secondary |
| Gate 2d Konovalenko routes (evaluation weeks, check-in 35%) | 1 | Observed, secondary |

In the σ_U > 0 cells, S2 and S3 use the carrier's mean times (the dispatcher's information), as L1 and L2 do.

## Hypotheses and decision rule (thresholds set now)

| ID | Statement | Threshold | Set |
|---|---|---|---|
| **H7a (primary)** | S1 recovers a material share of behind-schedule's shortfall | median G(S1) over the 14 cells at f = 0.35 (7 regimes × 2 noise levels), K/N = 5%, ≥ 0.25 | Gate 2c |
| **H7b (primary, no harm)** | S1 never does materially worse than L1 | eff(S1) ≥ eff(L1) − 0.02 in all 28 cells at K/N = 5% | Gate 2c |
| H7-rep | H7a replicates | median G(S1) over the 14 W cells at f = 0.35 ≥ 0.25 | Gate 3 W |
| H7-obs | The rule helps on observed outcomes | [V(W·B1) − V(B1)] / V(oracle), late-promise count, K/N = 5%: 95% bootstrap interval above 0 (2,000 resamples of routes, seed 2032) | Gate 4 LaDe; Gate 2d reported with the same statistic, no threshold |

| H7a | H7b | Outcome |
|---|---|---|
| holds | holds | **SUPPORT.** "Delay × remaining promises" is recommended as the model-free upgrade of behind-schedule triage when promise positions are unknown. It enters the ladder as rung L1w and the managerial implications. |
| holds | fails | **SUPPORT-WITH-CAVEAT.** As above, but the cells where it hurts are named and the recommendation is restricted. |
| fails | either | **NOT SUPPORTED.** The rule appears only as a theoretical remark. The shortfall of behind-schedule triage is attributed to heterogeneous slack, not to unequal sensitivities. |

H7-rep and H7-obs decide wording only: "replicates on a second sample" and "also on observed outcomes".

## Reported without thresholds

- Efficiencies of S1, S2 and S3 at 2%, 5% and 10% in every cell, next to L1, L2 and L3.
- G for each score, by regime.
- The share of the L2 (projection) advantage over L1 that each score recovers.
- In window regimes: the share of S1-selected routes with zero remaining promises, against that share for L1.
- Observed data: W·B2 as well as W·B1, and the minutes outcome.

## Verification before the run

1. For 50 random route-cells per simulation set, W_r recomputed from the route data with `promise_arrays` equals `U` (non-TW) or `windows_remaining` (TW) in the result file.
2. The script's own computation of eff(L1) and eff(L2) reproduces the published cell files (`gate2c_cells.csv`, `gate3_cells_A.csv`, `gate3_cells_B.csv`) to 1e-9.
3. Unit checks:
   - S2 / D = W when every c_ri = 0.
   - S3 = S1 when φ = 0.
   - λ ≥ 1 everywhere.
   - In the observed sets, the oracle and B1 values reproduce the published Gate 4 and Gate 2d efficiencies.

## Change log

| Date | Change | Reason |
|---|---|---|
| 2026-11-26 | Initial pre-registration | — |
| 2026-11-26 | Verification check 4 compared the Konovalenko B1 efficiency with the count outcome, but `gate2d_efficiency.csv` holds the late-minutes outcome only. The check was corrected to compare minutes with minutes (then 0 deviation), before the main run. No design change. | Bug in the verification script |
