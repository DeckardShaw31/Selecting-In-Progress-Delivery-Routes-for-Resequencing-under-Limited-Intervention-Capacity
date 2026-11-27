# Gate 4 report: second observed-outcome check on LaDe pickup tasks

Run 2026-11-26. The design is in `analysis/gate4/preregistration.md`, written after the files were downloaded and before any was opened. One deviation is logged there: the secondary score F-nowait was added at the schema check, before any score was computed.

- **Data:** LaDe-P (Cainiao), five Chinese cities, May–October. CC BY-NC 4.0 per the LaDe paper.
- **Code:** `code/resequencing/lade_io.py` (loading and cleaning) and `code/resequencing/gate4.py` (scores, metric, bootstrap).
- **Outputs:**
  - `cleaning.json`
  - `gate4_routes.csv`
  - `gate4_efficiency.csv`
  - `gate4_summary.json`
  - `gate4_within_cityday.csv`
  - `gate4_terciles.csv`
  - `gate4_descriptives.json`

## Decision

**Outcome under the pre-registered rule: SUPPORT.**

| Criterion (evaluation dates, check-in at 35% of tasks, top 5% of courier-days) | Value | Threshold | Result |
|---|---|---|---|
| D1 ratio V(behind schedule) / V(projection), late-task count | **0.777** | ≤ 0.85 | met |
| 95% bootstrap CI of [V(projection) − V(behind schedule)] / V(oracle), count | **[0.163, 0.195]** | above 0 | **met** |
| D2 (secondary): same test, late minutes | 0.814; CI [0.132, 0.164] | same | SUPPORT |

The behind-schedule gap now appears in observed outcomes on two independent datasets: delivery routes in Norway and Denmark (Gate 2d), and pickup tasks in China (this gate). Unlike in Gate 2d, the minutes outcome is also conclusive here. These data have far fewer timestamp errors, and the cleaning rule was set in advance.

## Data facts

- **Cleaning** (rules set in advance; whole courier-days excluded):
  - 751 courier-days were excluded for a task late by more than 24 h (1,089 tasks).
  - No task had a missing time, a window end before its start, or a pickup before acceptance.
- **Eligible:** 176,993 courier-days with at least 10 tasks (5.70 million tasks).
  - Calibration (first half of each city's dates): 85,260.
  - Evaluation: 91,733, of which 91,594 have at least one known remaining task at the 35% check-in.
- **Promises:** every window is exactly 120 min wide (5th and 95th percentiles both 120 min).
  - 1.7% of evaluation tasks are picked up after the window closes, and 26% before it opens.
  - 14.6% of courier-days have at least one known remaining task that ends up late (35% check-in).
- **Dynamic arrival:** at the 35% check-in, the median courier-day still has 30% of its tasks unaccepted, so they are unknown to the dispatcher. The median known remaining set is 10 tasks.
- **Fleet pace** (calibration, per city):
  - intercepts of 11.4–16.1 min per pickup and slopes of 1.1–7.5 min/km;
  - **R² 0.007–0.018**, so distance barely predicts the time between pickups.

## Results

Efficiency = share of the oracle's captured future lateness for the top-k courier-days ranked by each score. Selection is pooled across the five cities.

| Outcome | Check-in | k | Behind schedule (B1) | Currently late (B2) | Projection (F) | F-nowait | F-direct | F + route learning | Random |
|---|---|---|---|---|---|---|---|---|---|
| Late count | 0.35 | 5% | 0.600 | 0.401 | **0.772** | 0.775 | 0.792 | 0.839 | 0.069 |
| Late count | 0.35 | 10% | 0.526 | 0.362 | 0.766 | 0.769 | 0.710 | 0.820 | 0.114 |
| Late minutes | 0.35 | 5% | 0.612 | 0.407 | 0.752 | 0.753 | 0.808 | 0.800 | 0.060 |
| Late minutes | 0.35 | 10% | 0.551 | 0.377 | 0.739 | 0.741 | 0.736 | 0.778 | 0.103 |
| Late count | 0.50 | 5% | 0.703 | 0.419 | 0.851 | 0.852 | 0.839 | 0.861 | 0.068 |
| Late count | 0.50 | 10% | 0.664 | 0.374 | 0.846 | 0.849 | 0.822 | 0.866 | 0.110 |
| Late minutes | 0.50 | 5% | 0.680 | 0.415 | 0.854 | 0.856 | 0.899 | 0.875 | 0.059 |
| Late minutes | 0.50 | 10% | 0.647 | 0.386 | 0.857 | 0.859 | 0.864 | 0.872 | 0.101 |

**By city** (k = 5%, late count, check-in 0.35). Ratio = V(behind schedule) / V(projection):

| City | Courier-days | Behind schedule | Projection | Ratio |
|---|---|---|---|---|
| Chongqing | 19,899 | 0.662 | 0.853 | 0.776 |
| Hangzhou | 27,929 | 0.654 | 0.877 | 0.746 |
| Jilin | 4,929 | 0.169 | 0.545 | 0.310 |
| Shanghai | 19,411 | 0.676 | 0.828 | 0.817 |
| Yantai | 19,426 | 0.230 | 0.467 | 0.492 |

**The projection beats behind schedule in all 20 city cells** (5 cities × 2 outcomes × 2 check-ins). The ratios range from 0.13 to 0.85. The gap is largest in Jilin and Yantai, where behind schedule is barely better than random.

**Within city-day selection** (each city and date selects 5% of its courier-days, at least 1; check-in 0.35):
- late count: behind schedule 0.616 vs. projection 0.837;
- late minutes: 0.635 vs. 0.803.

The gap persists when selection is restricted to couriers working in the same city on the same day.

## Exploratory analyses (pre-registered without thresholds)

- **Window-width terciles** are uninformative: every window is 120 min wide, so the terciles split arbitrarily. Their ratios are 0.79, 0.76 and 0.79.
- **Share of the day's tasks still unknown at check-in.** The gap grows with this share. The D1 ratio is 0.88 in the lowest tercile (median 21% unknown), 0.74 in the middle (30%) and **0.54 in the highest (39%)**. When much of the day's work is still to arrive, lateness so far says little about the lateness to come, while the known remaining windows still do.

## What this adds to the paper

1. **A second observed dataset, independent of the first:** another country, another service type (pickups), another data provider. The direction and size of the gap match Gate 2d. Behind schedule reaches 73–83% of the projection there (clean outcomes) and 69–83% here (pooled cells).
2. **Prediction error on observed data.** The leg-level predictor is almost uninformative (R² ≈ 0.01), yet the projection still wins clearly. What matters is where the promises sit relative to the current time, not travel-time accuracy. This is the observed counterpart of the information-ladder claim: a promise-aware projection (L2) is enough where windows carry the exposure. LaDe is a pure window regime (ω = 1).
3. **Within-day selection** keeps the gap (count: 0.616 vs. 0.837).

## Limitations

- **Risk, not resequencing effect, is the target.** Gate 1 justified this on delivery routes; for pickups it is an assumption.
- **The courier's actual order is used to define "served".** The check-in is a share of the day's tasks, a definition that uses hindsight. Tasks accepted after the check-in are excluded from both the scores and the outcome.
- **Windows are uniform** (120 min), so the window-width dose analogue cannot be tested on LaDe.
- **Lateness is rare** (1.7% of tasks), so the top-5% selection concentrates on the few courier-days with future late tasks. Efficiencies are shares of a small total.

## Change log

| Date | Change |
|---|---|
| 2026-11-26 | Initial report |
