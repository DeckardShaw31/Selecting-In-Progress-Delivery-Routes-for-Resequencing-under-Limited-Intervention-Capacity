# Gate 2d report: behind-schedule vs. promise-aware triage on observed outcomes

Run 2026-11-26. Design: `analysis/gate2d/preregistration.md`, written before the data were downloaded; one logged post-hoc addition. Data: Konovalenko, Hvattum & Iversen (2026), Mendeley Data doi 10.17632/kkwgfvmtxn.1, `routes_performance.xlsx`, **CC BY 4.0**. 249,231 stop visits on 19,647 routes, 400 drivers, 32 weeks, Norway and Denmark. Code: `code/resequencing/gate2d.py` (pre-registered analysis), `gate2d_robustness.py` (follow-ups). Outputs: `gate2d_routes.csv`, `gate2d_efficiency.csv`, `gate2d_terciles.csv`, `gate2d_followups.csv`, `gate2d_summary.json`, `gate2d_followups_ci.json`.

## Decision

**Outcome under the pre-registered rule: CONSISTENT-WEAK.**

| Criterion (evaluation weeks, check-in at 35% of stops, top 5% of routes) | Value | Threshold | Result |
|---|---|---|---|
| D1 ratio V(behind schedule) / V(projection), future late minutes | **0.772** | ≤ 0.85 | met |
| 95% bootstrap CI of [V(projection) − V(behind schedule)] / V(oracle) | **[−0.020, 0.366]** | above 0 | **not met** (touches 0) |

## Data facts (descriptives, pre-registered)

- **Eligible routes:** 11,626 (at least 10 non-depot stops, complete times). Calibration weeks 0–15: 5,304 routes. Evaluation weeks 16–31: 6,322 routes.
- **Lateness:** 11% of stop visits arrive after the latest allowed time; 53% of eligible routes have at least one late stop. Median window width 240 min (5th–95th percentile 90–600).
- **Drivers deviate heavily:** 48% of non-depot stops are visited out of planned position, so a planned-order projection is handicapped.
- **Fleet pace fit (calibration):** 19.3 min per leg + 1.67 min per km; R² = 0.30 on 77,008 legs.
- **Data errors:** the maximum arrival time is ≈ 12.8 million minutes; 653 stop visits are "late" by more than 24 h; 7.7% of routes have arrival times that are not monotone in visit order. In the evaluation set the top 1% of routes carry 73% of all future late minutes, which is why the minutes-based interval is wide.

## Results

Efficiency = share of the oracle's captured future lateness (top-k routes ranked by the score).

**Primary outcome (future late minutes, all evaluation routes, pre-registered):**

| Check-in | k | Behind schedule (B1) | Currently late (B2) | Projection (F) | Projection + route learning | Random |
|---|---|---|---|---|---|---|
| 0.35 | 5% | 0.450 | 0.402 | **0.583** | 0.587 | 0.059 |
| 0.35 | 10% | 0.563 | 0.411 | **0.613** | 0.620 | 0.110 |
| 0.50 | 5% | 0.411 | 0.407 | **0.605** | 0.598 | 0.058 |
| 0.50 | 10% | 0.556 | 0.514 | **0.622** | 0.624 | 0.109 |

**Pre-registered secondary outcome (future late-stop count, all evaluation routes):** behind schedule 0.61–0.64; projection **0.74–0.81**; ratio 0.75–0.83 in all four cells. Bootstrap CI of the gap at check-in 0.35 and k = 5%: **[0.073, 0.174]**, above 0.

**Exploratory robustness (post hoc; excludes 66 routes with any stop late by > 24 h):** minutes outcome: behind schedule 0.49–0.55 vs. projection 0.63–0.75 (ratio 0.73–0.78); CI **[0.057, 0.237]**, above 0. Count outcome: unchanged from the secondary result.

**Exploratory window-width terciles** (minutes outcome, artifact-laden): no gap in the narrow and middle terciles, gap in the wide tercile (ratio 0.74). This is the opposite of the pre-stated expectation. With the minutes outcome dominated by a few erroneous routes, these terciles are not interpretable and are not used.

## Interpretation

- **On real routes with observed arrivals and real windows, the direction predicted by the simulations holds everywhere.** Projecting the remaining promises (planned order, fleet-average pace, no model) beats behind-schedule triage in all 12 cells (3 outcome definitions × 2 check-ins × 2 capacities). It captures about 1.1–1.5 times as much future lateness at the same capacity.
- **The magnitude is smaller than in the simulations.** Behind-schedule triage reaches about 73–83% of the projection's value on the clean outcomes (68–92% on the primary minutes outcome), against about 0.60 of risk triage in Gate 2c for observed windows. Possible reasons: much heavier driver deviation (48% out of order), very different routes (11 stops, mixed pickups and deliveries, depot returns), and different baselines (a crude projection here vs. a full posterior forecast in the simulations).
- **Route-specific learning adds almost nothing** (≤ 0.02), consistent with Gate 2's finding that most of the value comes from promise positions and variability, not from learning.
- **What this lets the paper say:** "An observed-outcome check on 6,322 real routes with actual arrivals and time windows reproduces the direction of the simulated gap. On the pre-registered secondary outcome (late-stop count), the gap is statistically clear [0.07, 0.17]; on the primary minutes outcome, data errors make the interval inconclusive." This is honest, and it materially strengthens the paper beyond "simulation calibrated to public data".

## Limitations

- **The target is risk, not recoverable lateness.** Realized future lateness is used because counterfactual resequencing is not computable here (no full distance matrix); Gate 1 justifies the substitution only under the simulation model.
- **The projection is crude.** It uses planned-order distances that may run from an already-served stop, and a fleet pace with R² of 0.30.
- **Behind schedule is measured against windows, not a plan.** The data have no planned arrival times, so the backward signal is lateness already incurred against windows, not delay against a planned schedule.
- **The primary outcome is contaminated by data errors in the source.**

## Change log

| Date | Change |
|---|---|
| 2026-11-26 | Initial report |
