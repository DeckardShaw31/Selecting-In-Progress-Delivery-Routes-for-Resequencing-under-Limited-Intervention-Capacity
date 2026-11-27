# Gate 2c pre-registration: confirmatory test of option 1 on independent routes

Written 2026-11-25, **before the Amazon training split was downloaded**. No training-split route has been seen. Gates 1 and 2 (evaluation split) motivate the hypotheses; this gate tests them on fresh routes with fixed thresholds. Changes after seeing results must be logged at the bottom.

## Hypotheses

- **C1 (replication of the gap).** Under day-definite + observed customer windows (MIXB), triage by delay so far captures at most 70% of the lateness reduction of predicted-risk triage (f = 0.35, base noise, K/N = 5%).
- **C2 (dose–response).** As the share of stops carrying exogenous customer-window promises rises, delay-triage efficiency falls monotonically.
- **C3 (mechanism).** Delay-triage efficiency is governed by **exposure endogeneity**: the more a route's exposure is created by its own delay (high corr(b, D) across routes), the better behind-schedule triage works. This measure was found post hoc in Gate 2 (Spearman 0.73 across 16 cells) and is now tested on independent routes.

## Data

- Amazon Last Mile Routing Research Challenge **training split** (`almrrc2021-data-training/model_build_inputs/`, 6,112 routes per the registry). Preprocessed exactly as the evaluation split (`code/resequencing/amazon_io.py`).
- **Sample:** 3,000 routes with at least 10 stops, drawn at random stratified by station (seed 2027), fixed before simulation. If fewer than 3,000 eligible routes exist, all are used.

## Execution model

Identical to Gates 1–2 (lognormal travel and service times with a route-day travel shock; posterior update from completed legs; SAA resequencing oracle with 48 in-sample and 256 out-of-sample scenarios; common random numbers). Noise: base and high. Check-in at f ∈ {0.35, 0.50}. Day-definite deadline = planned route end (with waiting) + 60 min.

## Promise regimes

Dose regimes (day-definite promise for stops without a window; customer window for stops with one):

| Regime | Stops with a customer window |
|---|---|
| DAY | none (0%) |
| MIXB | observed Amazon windows (≈ 6.6% on the evaluation split) |
| MIXB20 | observed + synthetic, total ≈ 20% of stops |
| MIXB50 | observed + synthetic, total ≈ 50% of stops |

Synthetic windows go to randomly chosen stops without a window (route-specific random stream, the same across noise levels and check-ins). Each gets an (opening offset, closing offset) pair relative to its planned no-wait ETA, drawn jointly from the empirical distribution of the 28,354 observed windows in the evaluation split (`analysis/gate2c/window_offsets_eval.npz`: median slack 4.2 h, 5th percentile 0.67 h). Waiting rules and the planned schedule are recomputed with the new windows.

Reference regimes (as in Gate 2, for the mechanism test): ETA, MIXA, TW.

## Measures (as in Gate 2)

- Information ladder L0–L4; value = sum of out-of-sample recoverable lateness over selected routes; efficiency = value / value of L4 risk triage.
- Exposure mass b_r = expected number of late promises (posterior, current sequence); exposure endogeneity = Pearson corr(b, D) across routes in a cell.
- Secondary (reported, no threshold): the share of expected late promises that are customer-window promises.

## Decision rule (thresholds set now)

- **C1:** MIXB, f = 0.35, base noise, K/N = 5%: V_L1/V_L4 ≤ 0.70.
- **C2:** in the dose series DAY → MIXB → MIXB20 → MIXB50, V_L1/V_L4 at K/N = 5% is non-increasing at every step in at least 3 of the 4 (f, noise) combinations, **and** V_L1/V_L4(DAY) − V_L1/V_L4(MIXB50) ≥ 0.20 at f = 0.35, base noise.
- **C3:** across all 28 cells (7 regimes × 2 check-ins × 2 noise levels), Spearman(V_L1/V_L4 at 5%, corr(b, D)) ≥ 0.70.

| Outcome | Action |
|---|---|
| C1, C2, C3 pass | **GO**: option 1 as a TR-E paper, with exposure endogeneity as the mechanism and the dose–response as the regime map. |
| C1 and C2 pass, C3 fails | **GO-NARROW**: the dose–response regime map stands, but the mechanism statistic is weak; the theory must rest on the identical-exposure proposition only. |
| C2 fails | **STOP option 1 as a TR-E paper**: report the comparison as an empirical note. |
| C1 fails, C2 passes | Gap smaller on fresh routes; GO-NARROW with the regime map as the main result. |

## Change log

| Date | Change | Reason |
|---|---|---|
| 2026-11-25 | Initial pre-registration | — |
