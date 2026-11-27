# Gate 4 pre-registration: second observed-outcome check on LaDe pickup tasks (E9)

Written 2026-11-26, **after the author downloaded the files and before any of them was opened**. Only file names and sizes have been checked:

| File | Size |
|---|---|
| `pickup_cq.csv` | 190.5 MB |
| `pickup_hz.csv` | 335.8 MB |
| `pickup_jl.csv` | 43.8 MB |
| `pickup_sh.csv` | 228.3 MB |
| `pickup_yt.csv` | 190.2 MB |

All five are in `data/raw/lade/pickup/`. What is known is the field list published with the dataset (Wu et al. 2024, LaDe; Hugging Face `Cainiao-AI/LaDe`). Any change made after seeing data must be logged at the bottom.

## Question

Gate 2d found, on observed Norwegian and Danish delivery routes, that a promise-aware projection beats behind-schedule triage. The gap was clear on the late-stop count and inconclusive on late minutes, because of timestamp errors.

On a second, independent dataset from another country and service type, does ranking couriers by lateness already incurred miss the courier-days that go on to incur the most promised-time lateness, relative to a forward-looking projection of the known remaining promises?

**Prediction from the theory.** Every LaDe pickup task carries a customer time window, so exposure is entirely customer-set (ω = 1, the TW regime). This is where the simulation shows the largest behind-schedule loss (0.41–0.54).

**Why Konovalenko v2 is not part of this gate.** A Konovalenko re-analysis with a count outcome and a cleaning rule was planned as part of E9. It is not included here because those numbers are already known from the Gate 2d follow-ups. Rerunning it would not be a fresh test.

## Data

**Source.** LaDe-P (pickup), five cities: Chongqing, Hangzhou, Jilin, Shanghai, Yantai. Licence: CC BY-NC 4.0 as stated in the LaDe paper (the Hugging Face metadata says apache-2.0; the paper is followed). Research use only; no redistribution.

**Fields used** (published names):
- courier ID and `ds` (date);
- `accept_time` (when the courier accepted the task);
- `time_window_start` and `time_window_end` (customer window);
- `pickup_time` (actual pickup);
- `lat`/`lng` (task location).

If a column name differs, it is mapped by meaning, and the mapping is logged.

**Parsing and distances.**
- Times are parsed as month-day clock times. A year boundary, if present, is resolved by the order of `ds`.
- Distances are haversine kilometres if coordinates lie within ±90°/±180°. Otherwise they are Euclidean in the given units, and the pace slope absorbs the unit.

**Route.** A route is one courier-day (city, courier, `ds`). Its tasks are ordered by actual pickup time; ties are broken by task row order.

**Cleaning rule (set now; applied to whole courier-days).** A courier-day is excluded if any task has:
- a missing or unparseable time or coordinate;
- `time_window_end` < `time_window_start`;
- `pickup_time` < `accept_time`;
- lateness (`pickup_time` − `time_window_end`) above 24 h.

**Eligibility.** At least 10 tasks after cleaning.

**Split by date, per city.** Each city's distinct `ds` values are sorted.
- The first half form the **calibration** set, used for pace parameters only.
- The second half form the **evaluation** set. All tests use evaluation courier-days only.

## Check-in and information

**When.** After m = round(f × n) tasks of the courier-day have been picked up, with f ∈ {0.35, 0.50}. The check-in time t0 is the pickup time of the m-th task, and the current position is that task's location.

**Known remaining tasks.** Tasks accepted at or before t0 and not among the first m picked up. Tasks accepted after t0 are unknown to the dispatcher at check-in and are excluded from both the scores and the outcome.

**Minimum.** A courier-day enters a check-in only if it has at least one known remaining task.

## Outcome (triage target)

- **Primary:** the count of known remaining tasks picked up after `time_window_end`. Following Gate 2d, the count is robust to timestamp outliers.
- **Secondary:** realized future late minutes, Σ (`pickup_time` − `time_window_end`)⁺ over the known remaining tasks.

As in Gate 2d, realized future lateness (risk) stands in for the effect of resequencing, justified by Gate 1. Gate 1 used delivery routes, so this is an assumption for pickups.

## Scores available at check-in

| Score | Definition | Role |
|---|---|---|
| **B1 late so far** (primary backward-looking) | Σ over the m picked-up tasks of (pickup − window end)⁺, minutes | Lateness already incurred (deviation-based exception signal) |
| B1c late count so far | Number of the m picked-up tasks that were late | Count version of B1 |
| B2 currently late | (pickup − window end)⁺ of the m-th task | "Behind now" |
| **F projection** (primary forward-looking) | Project the known remaining tasks from the current position at t0 in order of window end (ties: nearer first). Each leg takes a + b × distance; wait until window start if early. Score = Σ (projected pickup − window end)⁺ | Promise-aware projection with fleet pace (≈ L2); the order is one a dispatcher could compute |
| F-direct | Σ (t0 + a + b × distance from current position − window end)⁺ over known remaining tasks | Order-free lower bound on projected lateness |
| F-learn | As F, with the intercept replaced by a courier-day value from completed legs, shrunk toward a: a_r = (m·â_r + 10·a)/(m + 10) | Adds route-specific learning (≈ L4) |
| Oracle | The outcome itself | Upper bound |

**Pace.** (a, b) per city is the OLS fit of consecutive inter-pickup times on the distance between consecutive tasks. It uses all legs of calibration courier-days with 0 < inter-pickup time < the city's 99th percentile.

## Metric

- For a score S and capacity share k: V(S) = Σ of the outcome over the top-k share of evaluation courier-days ranked by S, pooled across the five cities. Ties are broken by (city, courier, `ds`) order.
- Efficiency = V(S) / V(Oracle).
- k ∈ {5%, 10%}.

## Decision rule (set now)

**D1 (primary):** f = 0.35, k = 5%, count outcome. Two conditions must both hold:
- V(B1) / V(F) ≤ 0.85;
- the 95% bootstrap interval of [V(F) − V(B1)] / V(Oracle) lies above 0. The interval uses 2,000 resamples of evaluation courier-days with seed 2030.

| Outcome | Interpretation for the paper |
|---|---|
| D1 passes | **SUPPORT**: the behind-schedule gap replicates on a second observed dataset (pickups, China). The observed-outcome evidence becomes two datasets in two service types. |
| V(B1)/V(F) ≤ 1 but D1 fails | **CONSISTENT-WEAK**: direction as predicted, not conclusive; report as such. |
| V(B1)/V(F) > 1 | **CONTRADICT**: on LaDe the backward signal does at least as well. Before submission, the simulated gap must be tempered, or explained by what differs (dynamic task arrival, pickup service). |

**D2 (secondary, same interval):** the same D1 test with the late-minutes outcome. Reported with its own SUPPORT / CONSISTENT-WEAK / CONTRADICT label. It does not change the D1 outcome.

## Reported without thresholds

- All score efficiencies for both outcomes, both check-ins and both k; random selection.
- Per-city D1 ratios.
- **Within city-day selection:** each (city, `ds`) selects max(1, round(0.05 N_d)) courier-days.
- **Exploratory dose analogue:** the D1 ratio within terciles of the median window width of the known remaining tasks.
- **Exploratory dynamic-arrival analysis:** the D1 ratio within terciles of the share of the day's tasks still unknown at check-in.
- **Descriptives:**
  - courier-days and tasks before and after cleaning, with the number excluded by each rule;
  - share of tasks late;
  - window widths;
  - share of tasks accepted after the check-in;
  - pace fit (a, b, R²) per city.

## Implementation

- **Scripts:** `code/resequencing/lade_io.py` (loading and cleaning) and `code/resequencing/gate4.py` (scores, metric, bootstrap), modelled on `gate2d.py`.
- **Order of work:**
  1. The schema check is the first step. Column names and time formats are printed, and any mapping is logged here.
  2. The decision is recorded before any exploratory analysis.

## Change log

| Date | Change | Reason |
|---|---|---|
| 2026-11-26 | Initial pre-registration | — |
| 2026-11-26 | **Schema check** (the first implementation step, done before any score or outcome was computed). All five files carry the published fields. The task ID is named `order_id` (published as package_id). Times are "MM-DD HH:MM:SS", and `ds` is MMDD. Coordinates are degrees, so haversine distances are used. Three sample rows show pickups before the window start. **Added a secondary score, F-nowait** (F without waiting for the window start). F stays primary as registered. | A projection that waits until window start may overstate projected lateness when couriers pick up early. Keeping F primary is the conservative choice, because it can only handicap the forward-looking score. |
