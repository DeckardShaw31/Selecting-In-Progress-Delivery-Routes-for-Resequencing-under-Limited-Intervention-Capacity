# Gate 2d pre-registration: observed-lateness check of behind-schedule vs. promise-aware triage

Written 2026-11-25, **before the dataset was downloaded**. Only the published field list is known (Konovalenko, Hvattum & Iversen 2026, *Data in Brief*, read in full). Changes after seeing data must be logged at the bottom.

## Question

Gates 1–2c are simulations. On real routes with **observed** arrival times and real time windows, does ranking in-progress routes by a backward-looking deviation signal miss the routes that go on to incur the most promised-time lateness, relative to a forward-looking projection of the remaining promises?

## Data

- Konovalenko, Hvattum & Iversen (2026), "Last-mile delivery route deviations dataset: planned vs. actual routes", Mendeley Data, doi 10.17632/kkwgfvmtxn.1, file `routes_performance.xlsx` (19.7 MB), licence **CC BY 4.0** (checked 2026-11-25 via the Mendeley public API).
- Fields used per stop: Route ID, Week ID, Depot flag, IndexP (planned position), IndexA (actual position), Arrived Time (actual arrival), Earliest Time, Latest Time (allowed arrival window), DistanceP, DistanceA (geodesic km to the previous stop in the planned/actual order). All times are minutes relative to a route-specific reference.
- **Eligible routes:** at least 10 non-depot stops, with arrival and latest times present for all non-depot stops.
- **Split by week:** weeks sorted by Week ID; routes in the first half of weeks form the **calibration** set (pace parameters only); routes in the second half form the **evaluation** set. All tests use the evaluation set only.

## Check-in

After m = round(f × n) non-depot stops have been served in actual order, f ∈ {0.35, 0.50}. Check-in time t0 = actual arrival time at the m-th served stop. Remaining stops are the unserved non-depot stops.

## Outcome (target of triage)

Realized future lateness of the route, Y = Σ over remaining stops of (Arrived − Latest)^+, in minutes. Secondary: the count of remaining stops arriving after Latest.

Because the dataset has no full distance matrix, the counterfactual effect of resequencing cannot be computed. Realized future lateness (risk) is used as the triage target. This is justified by Gate 1, where ranking by risk captured ≥ 96% of what ranking by recoverable lateness captured, but it is an assumption of this gate.

## Scores available at check-in

| Score | Definition | Role |
|---|---|---|
| **B1 late so far** (primary backward-looking) | Σ over served stops of (Arrived − Latest)^+ | Deviation-based exception signal (lateness already incurred) |
| B2 currently late (secondary backward) | (Arrived − Latest)^+ at the m-th served stop | "Behind schedule now" |
| **F projection** (primary forward-looking) | Project the remaining stops in **planned order** from t0: each leg takes a + b × DistanceP minutes (DistanceP = distance to the stop's planned predecessor); wait until Earliest if early; score = Σ (projected arrival − Latest)^+ | Promise-aware projection with fleet-average pace, no route-specific learning (≈ L2) |
| F-learn (secondary forward) | As F, but the per-leg intercept a is replaced by a route-specific value from the route's completed legs, shrunk toward a: a_r = (m·â_r + 10·a)/(m + 10) | Adds route-specific learning (≈ L4) |
| Oracle | Y itself | Upper bound |

Pace parameters (a, b) are the OLS fit of consecutive actual inter-arrival times on DistanceA over all completed legs of **calibration** routes, using legs with positive inter-arrival times below the 99th percentile. Known approximation: when a remaining stop's planned predecessor has already been served, DistanceP measures from that stop, not from the vehicle's position. Stated as a limitation.

## Metric

For a score S and capacity share k: V(S) = Σ of Y over the top-k share of evaluation routes ranked by S (ties broken by route order). Efficiency = V(S) / V(Oracle). k ∈ {5%, 10%}.

## Decision rule (set now)

- **D1 (primary):** f = 0.35, k = 5%: V(B1) / V(F) ≤ 0.85, **and** the 95% bootstrap interval (2,000 resamples of evaluation routes, seed 2028) of [V(F) − V(B1)] / V(Oracle) lies above 0.

| Outcome | Interpretation for the paper |
|---|---|
| D1 passes | **SUPPORT**: the behind-schedule gap appears in observed outcomes; the paper can report an observed-outcome check. |
| V(B1)/V(F) ≤ 1 but D1 fails | **CONSISTENT-WEAK**: direction as predicted, not conclusive; report as such. |
| V(B1)/V(F) > 1 | **CONTRADICT**: on real routes, the backward signal does at least as well; the simulated gap must be tempered or explained before submission. |

## Reported without thresholds

- All score efficiencies at both check-ins and both k; B2 and F-learn.
- **Exploratory dose analogue:** D1 ratio within terciles of the median width (Latest − Earliest) of the route's remaining windows. Expectation: a larger gap when windows are narrower.
- Descriptives: share of stops late, window widths, pace fit (a, b, R²), share of stops visited out of planned order.

## Change log

| Date | Change | Reason |
|---|---|---|
| 2026-11-25 | Initial pre-registration | — |
| 2026-11-26 | After the primary analysis (decision CONSISTENT-WEAK recorded first): (a) computed the pre-registered secondary outcome (late-stop count) with the same bootstrap; (b) **post-hoc, exploratory** robustness run excluding 66 evaluation routes with any stop late by > 24 h (`code/resequencing/gate2d_robustness.py`). Neither changes the recorded decision. | The source data contain timestamp errors (maximum arrival ≈ 12.8 million minutes; 653 visits late by > 24 h); the top 1% of routes carried 73% of all future late minutes, which made the minutes-based interval uninformative. |
