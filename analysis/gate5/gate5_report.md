# Gate 5 report: delay × delay sensitivity

Run 2026-11-26. The design is in `analysis/gate5/preregistration.md`, written before any Gate 5 score was computed. There is one logged fix to a verification check, made before the main run; there was no design change.

**Code:** `code/resequencing/gate5.py`. Verification passed:
- remaining promise weights match `promise_arrays` in 100 of 100 sampled route-cells;
- the published L1 and L2 efficiencies are reproduced exactly in 96 cells;
- unit checks pass;
- the published observed B1 efficiencies are reproduced exactly.

**Outputs** (in `analysis/gate5/`): `gate5_cells_*.csv`, `gate5_observed.csv`, `gate5_decision.json`.

**Run time:** 33 s. No re-simulation was needed.

## Decision

**Outcome under the pre-registered rule: NOT SUPPORTED.**

| Hypothesis | Value | Threshold | Result |
|---|---|---|---|
| **H7a:** median gap recovered by S1 (delay × remaining promises), Gate 2c, 14 cells at f = 0.35, K/N = 5% | **0.071** | ≥ 0.25 | **fails** |
| **H7b:** S1 never worse than L1 by more than 0.02 (28 cells) | worst −0.018 | ≥ −0.02 | holds |
| H7-rep: same as H7a on Gate 3 W | 0.046 | ≥ 0.25 | fails |
| H7-obs: LaDe, count outcome, 95% interval of [V(W·B1) − V(B1)]/V(oracle) | [0.000, 0.000] | above 0 | fails (uninformative; see below) |

**Pre-registered consequence.** The rule appears in the paper only as a theoretical remark. The shortfall of behind-schedule triage is attributed to heterogeneous slack (where the promises sit), not to unequal delay sensitivities.

## Results

**Gap recovered, by regime** (Gate 2c, f = 0.35, median over the two noise levels).

| Regime | ETA | MIXA | DAY | MIXB | MIXB20 | MIXB50 | TW |
|---|---|---|---|---|---|---|---|
| G(S1) | 0.11 | 0.08 | 0.08 | 0.01 | 0.04 | 0.06 | 0.20 |

- **TW is the only regime with a material gain** (0.18–0.23 across noise levels and check-ins). Most of it comes from the zero-promise routes. Behind schedule spends 5–16% of its selections on routes with no remaining customer window, which therefore carry no risk. Weighting by the remaining promise count gives those routes a score of zero.
- **Where a fifth to a half of the stops carry windows** (MIXB20 and MIXB50), weighting recovers 4–9% of the shortfall at base noise and the 35% check-in.
- **Against the projection,** in the window regimes S1 recovers a median of 7.5% of the projection's advantage over behind schedule. The projection, which knows where the promises sit, captures 0.88–0.97 of risk triage in these regimes at base noise.

**Secondary scores** (no threshold).

| Score | Median G (Gate 2c, f = 0.35) |
|---|---|
| S2, pace-extrapolated (perfect-learning sensitivity) | 0.12 |
| S3, model-based sensitivity (Proposition 7) | 0.08 |

Neither approaches the projection.

**Replication.**
- **Gate 3 W:** median G(S1) = 0.05. At base noise S1 is slightly worse than L1 in eight of 14 cells, by at most 0.03; the largest losses are in MIXB at both check-ins, 0.023 and 0.029. H7b was pre-registered for Gate 2c only; on this sample its 0.02 margin would be exceeded once.
- **High noise** gives larger gains (G up to 0.32 in MIXA at f = 0.50). There, delays are large and route-wide, and the remaining promise count matters more.

**Observed outcomes** (check-in at 35%, K/N = 5%).

| Data | Outcome | B1 | W·B1 | Projection F | 95% interval of the W·B1 gain |
|---|---|---|---|---|---|
| LaDe pickups | late count | 0.600 | 0.600 | 0.772 | [0.000, 0.000] |
| Konovalenko deliveries | late count | 0.613 | 0.622 | 0.736 | [−0.009, 0.035] |
| Konovalenko deliveries | late minutes | 0.450 | 0.557 | 0.583 | not computed (minutes outcome contaminated by timestamp errors, Gate 2d) |

- **LaDe is uninformative for this test.** Only 3.2% of courier-days are late at the check-in (2,911 courier-days, fewer than K = 4,580). Both scores therefore select every late courier-day plus the same zero-score ones, so the selections coincide in every bootstrap sample.
- **Konovalenko is informative** (14% of routes late), and the count gain is not distinguishable from zero.

## Interpretation

- **The theory's large-delay limit is correct, but its practical value is small.** In the regimes that matter, delays at the check-in are not large relative to the slack differences between routes, so the term $\tilde W_r(D_r-\tilde\mu_r)$ is dominated by the slack term $\tilde\mu_r$, not by the sensitivity.
- **This is consistent with Gates 2c and 3.** The mechanism statistic is ω, where the at-risk promises sit, and the projection, which uses slack positions, is what closes the gap.
- **The result sharpens the managerial message.** Counting remaining promises does not substitute for knowing where they sit. The one exception: under window-only promises, dropping routes with no remaining window is a cheap improvement.

## Limitations

- **Same simulation outputs.** The test reuses the Gate 2c and Gate 3 runs, so it inherits their simulator and promise regimes.
- **Approximate sensitivity in S3.** It uses the linear-Gaussian learning slope of Proposition 7, not the simulator's posterior. S2 uses the perfect-learning limit.
