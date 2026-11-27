# Pre-specified check of the average-case bound (Proposition 5 in `proof_audit.md`; Proposition 4 in the manuscript)

Written 2026-11-26. This is after the Proposition 5 derivation and the synthetic checks T5–T7 (all pass), and **before `code/resequencing/avg_bound.py` was run on any simulated cell**.

## What has and has not been seen

- **Gate 3 W cells (28):** $H_{\rm iso}$ was seen there and tracks the captured-risk loss at Spearman 0.97. That finding motivated Proposition 5, so these cells are **exploratory** for this check. The same holds for the Gate 3 σ_U > 0 cells (36), whose results were also inspected.
- **Gate 2c (28 cells) and Gate 3 Arm B (40 cells):**
  - The behind-schedule efficiencies of these cells are known.
  - $B_i$ is known on Gate 2c.
  - $H_{\rm iso}$ was computed in output files for both sets, but its values were not viewed.
  - No Proposition 5 quantity has been computed on them.
  - These two sets are the **pre-specified test sets**.
- **This is a weaker standard than the frozen-gate pre-registrations.** The data were generated before this note, and the motivating pattern came from closely related cells. The paper must describe it that way.

## Quantities (per cell; K = round(0.05 N); m = min(K, N − K))

- **Risk and delay:** risk ρ = `rho_min` (true posterior expected weighted lateness, L4*); delay D = `delay_plan_min` (L1).
- **Loss:** captured-risk loss 1 − V(S_D)/V(S_ρ), with the same stable tie-breaking as `analyze_gate3.py`.
- **B_quad_iso** (primary statistic): √(2m)‖ρ − ĝ(D)‖₂ / V(S_ρ), where ĝ is the least-squares isotonic fit.
- **B_tail_iso:** [T_m(e⁺) + T_m(e⁻)] / V(S_ρ) with the same ĝ.
- **B_tail_opt:** the tail bound minimized over nondecreasing ḡ (linear program).
- **B_inf:** 2m ε* / V(S_ρ), the Proposition 2′(i) term, for comparison.

## Predictions and thresholds

| ID | Prediction | Threshold | Sets |
|---|---|---|---|
| P0 | Every bound holds in every cell, and B_tail ≤ B_quad | all cells (a violation is a bug, to be fixed before reporting) | all four sets |
| **P1 (primary)** | B_quad_iso ranks cells by loss | Spearman ≥ 0.80 in **both** Gate 2c and Gate 3 Arm B | pre-specified sets |
| P2 | B_tail_iso and B_tail_opt rank cells by loss | Spearman ≥ 0.80 in both sets (each statistic separately) | pre-specified sets |
| P3 | Informativeness | no threshold: share of cells with B < 1, and median ratio of bound to loss, for each bound | all four sets |

## Consequences for the paper

| Result | Wording |
|---|---|
| P1 met | Proposition 5's quadratic term may be presented as a mechanism statistic alongside ω. The paper must say it was suggested by Gate 3 and checked on Gate 2c and Gate 3 Arm B. |
| P1 not met | Proposition 5 is presented as a sharper guarantee only (dominating Proposition 2′(i)); ω remains the mechanism statistic. |
| P0 fails | Stop; find the bug. |

The exploratory sets are reported alongside, labelled as exploratory.

## Results (run 2026-11-26 after this note was written; `analysis/theory/avg_bound/`)

| Statistic: Spearman with loss (share of cells with bound < 1; median bound/loss) | Gate 2c (pre-specified) | Gate 3 Arm B (pre-specified) | Gate 3 W (exploratory) | Gate 3 σ_U > 0 (exploratory) |
|---|---|---|---|---|
| **B_quad_iso (P1, primary)** | **0.65** (36%; 10.1) | 0.83 (35%; 7.4) | 0.76 (54%; 7.4) | 0.84 (44%; 5.6) |
| B_tail_iso (P2) | 0.94 (71%; 7.1) | 0.96 (68%; 4.5) | 0.95 (86%; 5.5) | 0.83 (89%; 3.6) |
| B_tail_opt (P2) | 0.96 (100%; 6.8) | 0.97 (95%; 4.2) | 0.97 (96%; 5.2) | 0.82 (97%; 3.4) |
| B_inf (Proposition 2′(i), comparison) | 0.29 (0%; 41) | 0.38 (0%; 30) | 0.33 (0%; 24) | 0.49 (0%; 25) |
| H_iso (statistic only) | 0.91 | 0.94 | 0.97 | 0.72 |

- **P0 met.** Every bound holds in every cell of all four sets (132 cells).
- **P1 not met.** The quadratic form reaches 0.65 on Gate 2c, below 0.80, and 0.83 on Arm B; the threshold needed both sets.
- **P2 met** for both tail statistics in both pre-specified sets.

**Decision under the rule above:** the tail-average bound is presented as a sharper guarantee, not as a mechanism statistic. The P2 results are reported as pre-specified secondary results, and ω remains the mechanism statistic.

**Why the quadratic form ranks less well:** the factor s_ρ / ρ̄_K differs between check-ins and noise levels, and the Cauchy–Schwarz step loses most when the misfit is spread out.
