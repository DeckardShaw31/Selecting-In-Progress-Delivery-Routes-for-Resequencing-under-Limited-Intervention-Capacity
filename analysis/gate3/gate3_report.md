# Gate 3 report: imperfect prediction, delay structure, within-day selection and the loss bound

Run 2026-11-26. The design is in `analysis/gate3/preregistration.md`, frozen before the main run, with no deviations.

The decision below was recorded from `analyze_gate3.py` output before any exploratory analysis. The figure code was edited afterwards for layout only, and H3's pre-registered consequence turns the delay-structure figure into a heat map. The edits changed no numbers.

**Code:**
- `code/resequencing/gate3.py` (simulation)
- `code/resequencing/gate3_verify.py` (four checks, all passed before the run)
- `code/resequencing/analyze_gate3.py` (decision rules, written before the run)

**Outputs** (in `analysis/gate3/`):
- `results/gate3_A.csv` (193,792 rows) and `results/gate3_B.csv` (121,120 rows)
- `gate3_cells_A.csv` and `gate3_cells_B.csv`
- `gate3_decision.json`
- `gate3_h4_replication_gate2c.csv`
- `gate3-prediction-error.png` and `gate3-delay-structure.png`

**Run:**
- Sample: 3,028 routes on 170 whole station-days, seed 2029.
- Arm A: 24,224 route-configurations in 2 h 05 min on 14 cores.
- Arm B: 30,280 route-configurations in 2 h 05 min.
- The background job reported exit code 1. That came from a trailing `tail -2` that Git Bash rejects; both arms logged "done", and every task key is present exactly once.

## Decision

**Outcome under the pre-registered rule: GO-NARROW (forecast).**

All primary tests use σ_U = 0.50, check-in f = 0.35 and K/N = 5%.

| Hypothesis | Value | Threshold | Result |
|---|---|---|---|
| H5a: V(L2) ≥ V(L1) in every cell with ω ≥ 0.50 | 5 of 5 qualifying cells | all | **holds** |
| H5b: V(L3)/V(L4*) ≥ 0.80 in all 12 cells | min 0.786 (MIXB, base); 0.792 (DAY, base) | all ≥ 0.80 | **fails** (2 of 12 cells below) |

Pre-registered consequences:
- **The projection claim stands.**
- **The practice rule changes.** The old rule, "a prior forecast is within about 12% of risk triage everywhere", is restated as conditional on prediction quality.

H5a in detail:

| Cell (σ_U = 0.50) | ω | V(L1) | V(L2) | L2 ≥ L1 |
|---|---|---|---|---|
| MIXB20, base | 0.580 | 0.712 | 0.812 | yes |
| MIXB50, base | 0.900 | 0.620 | 0.886 | yes |
| MIXB50, high | 0.748 | 0.797 | 0.872 | yes |
| TW, base | 1.000 | 0.448 | 0.929 | yes |
| TW, high | 1.000 | 0.449 | 0.884 | yes |

H5b in detail: V(L3)/V(L4*) per cell.

| Regime | ETA | DAY | MIXB | MIXB20 | MIXB50 | TW |
|---|---|---|---|---|---|---|
| base noise | 0.840 | **0.792** | **0.786** | 0.832 | 0.880 | 0.933 |
| high noise | 0.888 | 0.856 | 0.846 | 0.849 | 0.899 | 0.898 |

### Secondary hypotheses (wording only)

| Hypothesis | Value | Threshold | Result | Pre-registered consequence |
|---|---|---|---|---|
| H3: L1 efficiency rises with route-wide share φ | 4 of 4 series meet both conditions (gains +0.30 to +0.47; Spearman 1.00 in each) | ≥ 3 of 4 | **supported** | Delay structure is the second axis of the regime map; the heat map is a main exhibit |
| H3-D: L1 efficiency ≤ 0.80 at φ = 0.75, MIXB, base | 0.895 | ≤ 0.80 | **not met** | The loss at observed window shares is given as a range (φ = 0.75 to base value); the abstract leads with window-dominated settings |
| H4: Spearman(B_i, captured-risk loss of L1) over 28 W cells | 0.328 (0.335 against 1 − L1 efficiency) | ≥ 0.80 | **not met** | The bound is presented as a guarantee only; ω remains the empirical mechanism statistic |
| H6: within-day L1 efficiency, MIXB, f = 0.35, base | 0.661 | ≤ 0.80 | **met** | The gap is stated to persist within station-days |

### Bound checks (pre-registered verification 4)

- Proposition 2(i) holds in every cell of both arms.
- Proposition 2(ii) holds in every cell of both arms.
- The route-level Proposition 2(ii) inequality holds in every cell of both arms.
- In the Gate 2c replication (28 cells), Proposition 2(i) holds in all cells. Spearman(B_i, loss) = 0.286, which confirms that H4 fails there too.

## Results reported without thresholds

**W replication of Gate 2c.**
- On the new sample, behind-schedule efficiency in MIXB at f = 0.35 and base noise is **0.624**; Gate 2c gave 0.603.
- The other W cells reproduce the Gate 2c ordering: ETA and DAY ≥ 0.85; MIXB50 0.47 (base) and 0.76 (high); TW 0.32 (base) and 0.45 (high).

**σ_U curve** (Figure `gate3-prediction-error.png`; range over the 12 cells at each level).

| σ_U | Mean abs. error of leg means | L1 | L2 | L3 | L4f |
|---|---|---|---|---|---|
| 0 (W, f = 0.35) | 0 | 0.32–0.93 | 0.32–0.97 | 0.81–0.97 | 1.00 |
| 0.25 | about 20% | 0.36–0.85 | 0.41–0.94 | 0.82–0.96 | 0.94–0.99 |
| 0.50 | about 39% | 0.45–0.86 | 0.66–0.93 | 0.79–0.93 | 0.86–0.97 |
| 1.00 | about 77% | 0.62–0.87 | 0.78–0.88 | 0.77–0.88 | 0.76–0.87 |

- **The projection keeps its advantage.** V(L2) ≥ V(L1) in every cell with ω ≥ 0.50 at every σ_U: 10/10 at W, 5/5 at 0.25, 5/5 at 0.50 and 4/4 at 1.00.
- **Cost of misspecification, V(L4f)/V(L4*):** 0.94–0.99 at σ_U = 0.25, 0.86–0.97 at 0.50, and 0.76–0.87 at 1.00.

**Misspecified dispersion** (W, f = 0.35; dispatcher's σ scaled by 0.5 or 2).
- Halving the dispersion costs little: L4 is 0.99–1.02 of L4*, and L3 0.75–0.97.
- Doubling it costs more: L3 falls to 0.72–0.95 and L4 to 0.84–0.99. The lowest values are in MIXB and DAY.
- So a forecast that overstates uncertainty loses more than one that understates it.

**Delay structure** (Figure `gate3-delay-structure.png`; L1 efficiency, K/N = 5%).

| Regime | φ = 0.10 | 0.25 | 0.50 | 0.75 (D) | 0.90 |
|---|---|---|---|---|---|
| DAY, base | 0.81 | 0.88 | 0.91 | 0.93 | 0.94 |
| MIXB, base | 0.31 | 0.56 | 0.80 | 0.90 | 0.93 |
| MIXB50, base | 0.34 | 0.47 | 0.65 | 0.77 | 0.82 |
| TW, base | 0.27 | 0.36 | 0.47 | 0.54 | 0.57 |
| DAY, high | 0.92 | 0.91 | 0.92 | 0.93 | 0.94 |
| MIXB, high | 0.78 | 0.89 | 0.92 | 0.93 | 0.94 |
| MIXB50, high | 0.60 | 0.75 | 0.90 | 0.92 | 0.93 |
| TW, high | 0.40 | 0.52 | 0.62 | 0.66 | 0.71 |

- **ω falls as φ rises in mixed regimes** (MIXB, base: 0.82 at φ = 0.10 to 0.18 at 0.90). When delay is route-wide, the promises at risk shift from windows to the day-definite deadlines that delay already tracks.
- **The prior forecast stays at 0.85–0.98 in every Arm B cell.**

**Mechanism statistics** (28 W cells).

| Statistic | Spearman with captured-risk loss of L1 |
|---|---|
| B_i (pre-registered, H4) | 0.33 |
| H_iso (least-squares isotonic error) | 0.97 |
| ω | 0.92 |

- B_i ≥ 1 (vacuous) in 100% of cells, and B_ii ≥ 1 in 100%.
- Spearman(B_ii, captured-risk loss of L2) = 0.56.

**Within station-days** (W cells, K_d = max(1, round(0.10 N_d))).
- Behind-schedule efficiency: MIXB 0.66 (base) and 0.89 (high); MIXB20 0.55 and 0.84; MIXB50 0.50 and 0.76; TW 0.44 and 0.59.
- With K_d = 1 the numbers are similar (MIXB base 0.68).
- The prior forecast is 0.88–0.98 within days in every cell.
- Within-day L1 efficiency is 0.02–0.13 higher than pooled (largest in TW at base noise: 0.44 against 0.32). Pooling therefore overstates the gap modestly, but the gap persists.

## Exploratory observations (not pre-registered; for interpretation only)

1. **Behind schedule gains as prediction worsens; the forecasts lose.** In TW at base noise, L1 rises from 0.32 (W) to 0.62 (σ_U = 1.00); L3 and L4f fall to about 0.85.
   - Interpretation: with persistent errors in the carrier's means, delay so far is evidence that the route's true means exceed the carrier's. That evidence carries forward to future stops, and forecasts built on carrier means lack it.
   - The posterior forecast absorbs only the route-level part. The arc-level errors are independent by design.
   - This is a sharper, empirical form of the "learning" argument for H3.
2. **At σ_U = 1.00 the ladder collapses to a band of 0.62–0.88**, and the projection, prior forecast and posterior forecast are indistinguishable. Prediction quality is a real condition on the value of forecasting, though it never reverses the projection's advantage over behind schedule where windows dominate.
3. **The average-type isotonic error H_iso tracks loss almost perfectly (0.97).** The worst-case term ε* does not (0.33). The heterogeneity Proposition 2 names is the right quantity, but its L∞ version is driven by a few extreme route pairs.
   - This motivates a relative or average-case bound, an open Phase 6 item.
   - It is not a claim the paper can make as a tested mechanism.

## Consequences for the paper (applied in `manuscript/manuscript.tex` and `framework.md`)

- **Abstract and contributions:**
  - lead with window-dominated settings;
  - state the gap at observed window shares as a range: behind schedule forgoes 10–38% in MIXB at base total variance, from φ = 0.75 to the assumed delay structure;
  - present delay structure (φ) as the second axis of the regime map.
- **Practice rule:** "use the projection where windows dominate at-risk promises" stands under prediction errors up to about 77% per leg. "A prior forecast captures most of the value" becomes conditional on prediction quality: about 0.79–0.93 at a 39% leg-mean error, below 0.80 in two of 12 cells.
- **Proposition 2:** presented as a guarantee (it is never violated, and it is tight in the worst case) but not as a predictor of loss across cells; ω remains the mechanism statistic.
- **Within station-days:** the gap persists (H6 met).

## Limitations

- **Arc-level independence.** The prediction errors are independent across arcs by design, which is the hardest case for learning. Spatially correlated errors would favour the posterior forecast.
- **Oracle resequencing.** The dispatcher never resequences with the wrong model; value is measured with oracle resequencing.
- **The calibrated φ = 0.75 comes from one Norwegian–Danish fleet** (Konovalenko et al. 2026) and may not transfer.
