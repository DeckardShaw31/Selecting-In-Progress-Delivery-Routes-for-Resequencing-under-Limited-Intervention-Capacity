# Gate 2 report: behind schedule vs. at risk (option 1)

Run 2026-11-25. Design and thresholds: `analysis/gate2/preregistration.md` (written before the run; one implementation fix logged before the main run). Code: `code/resequencing/gate1.py` (regimes MIXA/MIXB and ladder scores added; Gate 1 outputs verified unchanged), `analyze_gate2.py`. Results: `analysis/gate2/results/gate2_main.csv` (48,832 route-configurations: 3,052 routes × 4 promise regimes × 2 check-ins × 2 noise levels; 38 min on 14 cores), cell table `gate2_cells.csv`.

## Decision

**Outcome under the pre-registered rule: REVISE.**

| Criterion | Value | Threshold | Result |
|---|---|---|---|
| G2a: delay-triage value / risk-triage value, MIXB (day-definite + customer windows), f = 0.35, base noise, K/N = 5% | **0.542** | ≤ 0.70 | **pass** |
| G2b(i): Spearman(delay-triage efficiency, corr(D, ρ)) across 16 cells | 0.591 | ≥ 0.70 | fail |
| G2b(ii): Spearman(corr(D, ρ), 1/√(1+CV(b)²)) across 16 cells | 0.453 | ≥ 0.70 | fail |

The practical gap is real: in the primary realistic regime, triaging by "behind schedule" removes only 54% of the lateness that triaging by predicted risk removes. But the pre-registered mechanism (heterogeneity of exposure mass across routes, summarized by CV(b)) does **not** explain where the gap appears. Without a confirmed mechanism, the pre-registered consequence is an empirical-only paper, which is weaker for TR-E, unless a re-specified mechanism is confirmed on independent data (see "Recommended next step").

## Results (K/N = 5%; values are lateness removed relative to L4 risk triage)

| Regime | Check-in | Noise | L0 random | **L1 delay so far** | L2 projection | L3 + variability | L4 risk | L4 / oracle |
|---|---|---|---|---|---|---|---|---|
| ETA | 0.35 | base | 0.05 | **0.92** | 0.28 | 0.90 | 1 | 0.99 |
| ETA | 0.35 | high | 0.06 | **0.88** | 0.59 | 0.89 | 1 | 0.99 |
| ETA | 0.50 | base | 0.05 | **0.96** | 0.68 | 0.95 | 1 | 0.99 |
| ETA | 0.50 | high | 0.06 | **0.89** | 0.81 | 0.93 | 1 | 0.99 |
| MIXA | 0.35 | base | 0.06 | **0.84** | 0.65 | 0.93 | 1 | 0.95 |
| MIXA | 0.35 | high | 0.06 | **0.86** | 0.64 | 0.87 | 1 | 0.99 |
| MIXA | 0.50 | base | 0.06 | **0.83** | 0.85 | 0.97 | 1 | 0.96 |
| MIXA | 0.50 | high | 0.06 | **0.88** | 0.82 | 0.93 | 1 | 0.99 |
| MIXB | 0.35 | base | 0.07 | **0.54** | 0.81 | 0.88 | 1 | 0.99 |
| MIXB | 0.35 | high | 0.07 | **0.81** | 0.69 | 0.84 | 1 | 0.99 |
| MIXB | 0.50 | base | 0.06 | **0.57** | 0.84 | 0.93 | 1 | 0.99 |
| MIXB | 0.50 | high | 0.07 | **0.86** | 0.80 | 0.89 | 1 | 0.98 |
| TW | 0.35 | base | 0.06 | **0.30** | 0.95 | 0.95 | 1 | 1.00 |
| TW | 0.35 | high | 0.07 | **0.39** | 0.89 | 0.93 | 1 | 1.00 |
| TW | 0.50 | base | 0.06 | **0.33** | 0.96 | 0.97 | 1 | 0.99 |
| TW | 0.50 | high | 0.06 | **0.46** | 0.93 | 0.95 | 1 | 0.99 |

L1 = rank by delay so far (realized minus planned completion at check-in, regime-consistent plan). L2 = lateness of the deterministic projection of the current sequence with mean times from the realized position and time (no forecasting model). L3 = expected lateness under the prior model (variability, no route-specific learning). L4 = expected lateness under the posterior (Gate 1 risk). Oracle = true recoverable lateness. At K/N = 2% and 10%, delay triage captures 0.67–0.70 and 0.54–0.56 in MIXB base, and 0.31–0.33 and 0.36–0.38 in TW base (`gate2_cells.csv`).

Figures: `gate2-information-ladder.png` (ladder by regime, f = 0.35, base) and `gate2-exposure-mechanism.png` (pre-registered mechanism test; points far off the predicted line).

## What the information ladder shows

1. **Behind-schedule triage is near-optimal only when promises are anchored to the plan.** With ETA-based promises it captures 0.88–0.96 of risk triage; with customer windows only 0.30–0.46; mixed regimes lie in between (MIXA 0.83–0.88; MIXB 0.54–0.86).
2. **Under customer windows, knowing where the promises are is almost all you need.** A deterministic projection with no forecasting model (L2) captures 0.89–0.96 under TW, closing 82–94% of the delay-to-risk gap (pre-stated rule: ≥ 80% → "promise-aware projection suffices").
3. **Under plan-anchored promises the projection fails and the delay itself is informative.** Under ETA promises L2 captures only 0.28–0.81, worse than delay so far: most routes look on time when projected with mean times, while the delay so far carries information about the route-level slowdown. A probabilistic forecast (L3) is needed to match or beat delay triage, and route-specific learning (L4) adds the last 5–11%.
4. **In the primary regime (MIXB, base noise) the answer is "both":** the projection closes 58–64% of the gap and the prior forecast 73–84%.
5. **Risk triage stays near the oracle** (L4 / oracle ≥ 0.95 in every cell), consistent with Gate 1.

## Why the pre-registered mechanism failed (exploratory, not part of the decision)

`analysis/gate2/exploratory_mechanism.csv`. CV(b) over all routes is large in every regime (2.9–5.9), because most routes carry no exposure; it does not separate the regimes (Spearman with delay-triage efficiency 0.12). The linearized model behind the prediction assumed exposure mass independent of delay. That assumption holds under customer windows but fails under ETA promises, where a route becomes exposed *because* it is late. A post-hoc measure of this endogeneity, corr(b, D) across routes, tracks delay-triage efficiency with Spearman 0.73 across the 16 cells; rank agreement between delay and risk (Spearman ρ–D) tracks it at 0.83–0.84. These are exploratory, computed on the same data, and must be confirmed on fresh routes before being claimed.

Candidate re-specified mechanism: **exogenous vs. endogenous exposure.** When promises are derived from the plan (ETA + buffer), exposure is generated by the delay itself, so "behind schedule" is a sufficient triage statistic up to variability. When promises are set by customers (windows), exposure is positional and exogenous to the delay, so delay-based triage misses most of the value, and a promise-aware projection recovers it.

## Recommended next step (for the author to decide)

A confirmatory gate on **independent data**: the Amazon **training split** (6,112 routes, 2.3 GB, not yet used), pre-registered before download, with:

- the re-specified mechanism (exposure endogeneity, e.g., corr(b, D) or a structural measure computed from promise positions alone) as the hypothesis, with thresholds set in advance;
- a **dose–response design**: vary the share of stops with exogenous promises (0%, the observed 6.6%, 20%, 50%, by assigning synthetic customer windows drawn from the observed window distribution) at fixed noise, and test whether delay-triage efficiency falls monotonically with that share. A monotone regime map is exactly the kind of result TR-E reviewers find convincing.

If confirmed, the paper becomes: *"Behind schedule is not at risk: promise structure decides what information capacity-limited last-mile interventions need"*. That framing has an analytical part (Proposition: delay ranking is risk-optimal when exposure profiles are identical; it degrades with exogenous positional exposure), a decomposition (positions vs. variability vs. learning), and a regime map. If not confirmed, the defensible output is an empirical regime comparison, better suited to a practice-oriented venue than TR-E.

## Limitations

Same as Gate 1: a perfectly specified predictor; assumed stochastic parameters; simulated execution on Amazon routes (no arrival timestamps); selection pooled across station-days; free complete resequencing. In addition, the MIXA/MIXB promise constructions are analyst-designed (logged in the pre-registration) and the day-definite deadline is planned end + 60 min.

## Change log

| Date | Change |
|---|---|
| 2026-11-25 | Initial report |
