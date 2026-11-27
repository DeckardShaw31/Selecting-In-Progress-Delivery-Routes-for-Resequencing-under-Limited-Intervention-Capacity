# Gate 2c report: confirmatory test of option 1 on independent routes

Run 2026-11-25. Design and thresholds: `analysis/gate2c/preregistration.md`, written before the training split was downloaded (no changes after the run). Data: Amazon training split, 6,112 routes, **zero route IDs shared with the evaluation split**; pre-registered stratified sample of 2,999 routes (seed 2027). Results: `analysis/gate2c/results/gate2c_main.csv` (83,972 route-configurations: 7 promise regimes × 2 check-ins × 2 noise levels; 45 min on 14 cores). Cell table: `gate2c_cells.csv`. Figure: `gate2c-dose-response-mechanism.png`.

## Decision

**Outcome under the pre-registered rule: GO-NARROW** (dose–response regime map confirmed; primary mechanism statistic weak).

| Criterion | Value | Threshold | Result |
|---|---|---|---|
| C1: delay triage / risk triage, MIXB (day-definite + observed windows), f = 0.35, base, K/N = 5% | **0.603** | ≤ 0.70 | **pass** (Gate 2: 0.542) |
| C2: dose series monotone (non-increasing) in ≥ 3 of 4 (f, noise) combinations, and DAY − MIXB50 ≥ 0.20 at f = 0.35 base | **4 / 4**; drop **0.494** | ≥ 3; ≥ 0.20 | **pass** |
| C3: Spearman(delay-triage efficiency, corr(b, D)) across 28 cells | 0.640 | ≥ 0.70 | fail |

Pre-registered consequence (GO-NARROW): the regime map stands as the main result; the theory must not lean on the exposure-endogeneity statistic. It rests on the identical-exposure proposition (delay ranking equals risk ranking when all routes share one exposure function), with the dose–response as the empirical regime map.

## Dose–response (delay-so-far triage relative to predicted-risk triage, K/N = 5%)

| Check-in | Noise | DAY (0% windows) | MIXB (observed, 7.1%) | MIXB20 (20%) | MIXB50 (50%) |
|---|---|---|---|---|---|
| 0.35 | base | 0.915 | 0.603 | 0.531 | 0.422 |
| 0.50 | base | 0.948 | 0.634 | 0.526 | 0.437 |
| 0.35 | high | 0.898 | 0.878 | 0.847 | 0.793 |
| 0.50 | high | 0.953 | 0.931 | 0.914 | 0.826 |

Behind-schedule triage is near-optimal when every promise is anchored to the plan (day-definite: 0.90–0.95). It loses 37–58% of the achievable lateness reduction once customer-window promises are present at observed or higher shares (base noise). When route-wide shocks dominate (high noise), the delay so far carries that shock and the loss is much smaller (17–21% at 50% windows): **two axes** govern the regime map: where the promises are (exogenous windows vs. plan-anchored) and what drives delay (route-wide vs. stop-level shocks).

## Secondary measure (pre-registered, reported without a threshold)

Share of a cell's expected late promises that are customer-window promises vs. delay-triage efficiency: **Spearman −0.921** across all 28 cells (−0.983 within base noise, −0.874 within high noise). For comparison, the raw share of stops with windows: −0.399; the primary endogeneity statistic: +0.640. The mechanism is therefore better described as **where the exposure sits** (on exogenous window promises vs. on plan-anchored promises) than as a correlation between exposure and delay. This matches the job-shop result of Baker (1984): the way due dates are assigned changes which priority rule works, and progress-against-milestone rules work when milestones reflect the plan's work content. Because it was pre-registered as secondary, it supports but does not change the decision.

## Information ladder on fresh routes (K/N = 5%, all 28 cells)

- **Promise-aware projection, no forecasting model (L2):** 0.88–0.99 wherever windows carry much of the exposure (MIXB, MIXB20, MIXB50, TW at base noise); but 0.16–0.56 under day-definite or ETA promises at base noise, where it is far worse than the delay itself.
- **Adding generic variability (L3):** ≥ 0.885 in every cell.
- **Route-specific learning (L4) vs. oracle recoverability:** 0.947–0.999 (Gate 1's near-optimality of risk ranking replicates).
- **Random selection:** 0.05–0.18.

## How the reading (all 23 closest papers on file) positions this result

- Forward-looking vs. backward-looking priority is established in job shops (Vepsäläinen & Morton 1987, read: FCFS called "backward-looking"; look-ahead slack rules win). The interaction between due-date setting and priority rules is also established (Baker 1984, read). Option 1 must not claim either.
- A supply-chain disruption system already displays backward- and forward-looking tour status (Meyer et al. 2018, read) but never quantifies which triage works or when; deviation-threshold alerting is the documented practice (Otto 2003, read).
- The recent dynamic-routing competitor (Liu, Ma & Ma 2026, read) re-plans all routes each epoch; no budget, no triage.
- **Defensible contribution:** the routing analogue of the due-date × priority-rule interaction for *capacity-limited interventions on in-progress routes*, where each route carries many promises at different positions and delays propagate along the remaining sequence: (i) a confirmed regime map (promise structure × delay-shock structure) on two independent samples of real routes; (ii) an information ladder showing when a promise-aware projection suffices and when a probabilistic forecast is needed; (iii) the identical-exposure proposition with the exposure-share statistic as its empirical counterpart.

## Limitations

- **Simulated outcomes.** The predictor is perfectly specified, the stochastic parameters are assumed, and lateness is simulated (the Amazon data have no arrival timestamps).
- **Synthetic window promises.** They are drawn from the empirical window distribution, but their placement on stops is random.
- **Selection pooled across station-days.**
- **The intervention is a free, complete reorder.**
- **An observed-lateness check is still missing.** The Konovalenko et al. (2026) dataset (actual arrivals and time windows on 19,497 real routes) could test the triage comparison against observed lateness; its data licence still needs checking.

## Change log

| Date | Change |
|---|---|
| 2026-11-25 | Initial report |
