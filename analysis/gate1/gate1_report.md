# Gate 1 report: recoverable vs. locked-in lateness on Amazon routes

Run 2026-11-25. Design and thresholds: `analysis/gate1/preregistration.md` (written before any simulation). Code: `code/resequencing/` (`sim.py`, `gate1.py`, `analyze_gate1.py`, tests in `test_sim.py`). Raw results: `analysis/gate1/results/`.

## Decision

**Outcome under the pre-registered rule: REVISE (risk ranking is near-optimal).**

| Criterion (f = 0.35, base noise) | ETA promises | TW promises | Threshold | Result |
|---|---|---|---|---|
| G1a: recoverable share R = ΣΔ / Σρ | 0.356 | 0.977 | ≥ 0.10 | **pass** (both) |
| G1b: V_risk / V_rec at K/N = 10% | 0.993 | 1.000 | ≤ 0.85 | **fail** (both) |

Resequencing the remaining stops removes a large share of expected lateness, so the intervention lever is material. But choosing routes by predicted lateness (risk) captures 99–100% of what choosing by recoverable lateness captures, so the core mechanism of the current framing (H1: "fix the recoverable routes, not the latest ones") is **not supported** on these routes. The pre-registered response is to recentre the paper. Options are in the last section.

## What was run

- **Data:** all 3,052 routes in the evaluation split of the 2021 Amazon Last Mile Routing Research Challenge (the registry page says 3,072; the file has 3,052). The routes average 142 stops (median 145, range 19–225) and last about 7.5 h nominally. Only 6.6% of stops have a customer time window, and the median window is 8 h wide; 92% of routes have at least one windowed stop.
- **Model (assumed):** lognormal travel times with a route-day shock and leg noise; lognormal service times. The dispatcher knows the model and updates the route shock from the legs completed so far (a perfectly specified predictor). This makes Gate 1 an upper-bound check, not a test under realistic prediction error.
- **Promises:** TW = Amazon customer windows (observed); ETA = planned arrival + 60 min (derived from the dispatched sequence).
- **Intervention:** free reordering of all remaining stops by SAA local search (48 scenarios), evaluated out of sample (256 scenarios, common random numbers).
- **Grid:** 2 promise regimes × 2 check-in points (35% and 50% of stops served) × 3 noise levels = 36,624 route-configurations, 571 s on 14 cores.
- **Verification:** the hand-computed example, a pure-Python reference evaluator and brute-force optima (local search optimal in 78 of 80 random instances) all pass.

## Pre-registered results

| regime | check-in | noise | late routes | mean ρ (min) | mean Δ (min) | R | Δ < 0 out of sample | Spearman(ρ, Δ) | V_risk/V_rec 5% | 10% | 20% | top-10% overlap |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ETA | 0.35 | low | 2.5% | 1.6 | 0.7 | 0.431 | 0.5% | 0.45 | 0.997 | 0.997 | 0.997 | 0.76 |
| ETA | 0.35 | base | 18.6% | 26.1 | 9.3 | 0.356 | 5.4% | 0.36 | 0.987 | 0.993 | 0.991 | 0.80 |
| ETA | 0.35 | high | 49.0% | 182.9 | 50.8 | 0.278 | 14.8% | 0.47 | 0.989 | 0.988 | 0.992 | 0.89 |
| ETA | 0.50 | low | 1.8% | 1.6 | 0.8 | 0.494 | 0.2% | 0.52 | 1.000 | 1.000 | 1.000 | 0.83 |
| ETA | 0.50 | base | 14.9% | 24.5 | 8.2 | 0.336 | 4.7% | 0.31 | 0.990 | 0.992 | 0.991 | 0.70 |
| ETA | 0.50 | high | 41.3% | 164.9 | 43.0 | 0.261 | 12.4% | 0.44 | 0.988 | 0.992 | 0.993 | 0.93 |
| TW | 0.35 | low | 13.2% | 4.3 | 4.2 | 0.979 | 0.1% | 0.95 | 0.996 | 0.999 | 1.000 | 0.99 |
| TW | 0.35 | base | 16.3% | 5.4 | 5.3 | 0.977 | 0.1% | 0.92 | 0.996 | 1.000 | 1.000 | 0.99 |
| TW | 0.35 | high | 23.0% | 9.8 | 9.5 | 0.968 | 0.7% | 0.88 | 0.998 | 0.999 | 0.999 | 0.99 |
| TW | 0.50 | low | 11.7% | 4.0 | 3.8 | 0.951 | 0.0% | 0.96 | 0.987 | 0.999 | 1.000 | 0.99 |
| TW | 0.50 | base | 14.0% | 5.1 | 4.8 | 0.944 | 0.1% | 0.92 | 0.989 | 0.999 | 1.000 | 0.98 |
| TW | 0.50 | high | 19.1% | 9.1 | 8.5 | 0.934 | 0.4% | 0.90 | 0.988 | 0.998 | 1.000 | 0.98 |

"Late routes" = expected lateness above 0.5 min. ρ = expected lateness of the current sequence; Δ = expected lateness removed by resequencing; both are out of sample and averaged over all routes (zeros included). Selection is pooled across the sample, because station-days in the evaluation split are incomplete (median 2 routes per station-day).

Figures: `gate1-rho-vs-delta_prereg.png` (route-level ρ against Δ) and `gate1-capacity-curve_prereg.png` (V_risk / V_rec against K/N).

## Why risk ranking is near-optimal (diagnosis)

1. **Lateness is concentrated.** Under ETA promises (base, f = 0.35), the 5% latest routes carry 96% of all expected lateness and 99% of all recoverable lateness; under TW promises, 79% and 79%. Any rule that finds the latest routes finds almost all the recoverable lateness.
2. **Recoverability rises with lateness, the opposite of the "locked-in" hypothesis.** Among late ETA routes, the recoverable share by lateness quintile is −0.24, 0.00, 0.12, 0.27, 0.37. Badly late routes have many late stops whose order can still be improved; mildly late routes have little to gain. Under TW promises almost all lateness is recoverable (Δ ≈ ρ), because the few windowed stops can simply be visited earlier. In both cases ρ and Δ are close to comonotone where it matters, which is exactly the condition of Proposition 1 under which risk ranking is optimal.
3. **The two top-K sets differ, but the differences are worth little.** At 10% capacity the sets overlap only 70–93% under ETA promises, yet the value differs by ≤1.2%: the routes swapped in and out have small Δ.
4. **Pricing interventions by disruption does not change it (quick check on the base cells).** With a budget on total stops moved (10% of all moves) instead of routes, ordering by risk still gets 96.5% (ETA) and 93.8% (TW) of the best ordering (by Δ per moved stop).

## Findings that survive and matter

- **The lever is material.** Resequencing at a mid-route check-in removes about one third of expected lateness against ETA promises (R = 0.26–0.49 across noise levels) and 93–98% of lateness against customer time windows.
- **Selection matters a lot, and prediction drives it.** Random selection captures about K/N of the achievable reduction (2–10%). Ranking by predicted risk captures 97–100%. A no-prediction "behind schedule" rule (rank by delay so far) captures 78–99% under ETA promises but only **4–44% under customer time windows**, because exposure depends on where the windowed stops sit in the remaining route, not on how late the driver already is (`policy_ratios_prereg.csv`).
- **Small capacity is enough.** Because lateness is concentrated, touching 5% of routes captures almost all recoverable lateness when routes are chosen by predicted risk.
- **Resequencing can backfire.** With 48 planning scenarios, the new sequence is worse out of sample on 5–15% of routes under ETA promises (more under high noise), mostly mildly late routes. This optimizer's-curse effect appears even with a perfectly specified model; real predictors will make it worse.
- **Full resequencing is disruptive.** On late routes the oracle moves a median of 25–33 stops (ETA) and 62–64 stops (TW). Drivers are unlikely to accept changes of this size, so realistic interventions are much more limited.

## Exploratory checks (not part of the decision rule)

Logged in the pre-registration change log before these runs; base noise, all 3,052 routes.

| Variant | regime | f | late routes | R | V_risk/V_rec 5% | 10% | 20% |
|---|---|---|---|---|---|---|---|
| 1-hour ETA windows (±30 min) | ETA | 0.35 | 53.1% | 0.344 | 0.987 | 0.990 | 0.993 |
| 1-hour ETA windows (±30 min) | ETA | 0.50 | 45.6% | 0.331 | 0.980 | 0.992 | 0.994 |
| Route-level service shock (sP = 0.20) | ETA | 0.35 | 36.5% | 0.300 | 0.964 | 0.994 | 0.998 |
| Route-level service shock (sP = 0.20) | ETA | 0.50 | 30.2% | 0.288 | 0.962 | 0.993 | 0.999 |
| Route-level service shock (sP = 0.20) | TW | 0.35 | 20.2% | 0.977 | 0.998 | 0.999 | 1.000 |
| Route-level service shock (sP = 0.20) | TW | 0.50 | 17.9% | 0.961 | 0.992 | 0.998 | 1.000 |
| 1-hour windows + service shock | ETA | 0.35 | 53.3% | 0.271 | 0.933 | 0.971 | 0.992 |
| 1-hour windows + service shock | ETA | 0.50 | 46.9% | 0.259 | **0.864** | 0.951 | 0.995 |
| Count KPI (number of late stops) | ETA | 0.35 | 11.9% | 0.273 | 0.986 | 0.983 | 0.978 |
| Count KPI (number of late stops) | ETA | 0.50 | 9.4% | 0.317 | 0.991 | 0.983 | 0.978 |
| Count KPI (number of late stops) | TW | 0.35 | 11.6% | 0.949 | 0.985 | 0.990 | 1.000 |
| Count KPI (number of late stops) | TW | 0.50 | 11.1% | 0.921 | 0.984 | 0.980 | 1.000 |

For the count KPI, "late routes" means more than 0.05 expected late stops, and R is in late stops rather than minutes.

- Tighter promises alone, or systematic service-time delays alone, raise lateness but leave the conclusion unchanged.
- **The one place where the rankings start to diverge is when lateness is widespread rather than concentrated** (1-hour windows plus a route-level service shock: about half of all routes late). There, risk ranking falls to 0.86–0.93 of recoverability ranking at 5% capacity, and to 0.95–0.97 at 10%. This is consistent with Proposition 1: comonotonicity holds at the top of a concentrated lateness distribution, and erodes as lateness spreads across many routes. It is a possible boundary condition (for example peak season or bad weather), but the effect is modest and still above the pre-registered 0.85 threshold at 10% capacity.
- On the on-time-rate (count) KPI, risk ranking still captures 98–99% at 5–10% capacity (0.90 at 2% under ETA promises), and the "behind schedule" rule again captures only 10–24% under customer time windows.

## Limitations of this gate

- The predictor is perfectly specified (the dispatcher knows the generating model), so H2 (noise and the reversal condition) is untested here.
- All stochastic parameters are assumed; the Amazon data have no arrival timestamps.
- Selection is pooled across routes from different station-days.
- The intervention is a free, complete reorder, with no limit on stops moved and full driver adoption.
- The objective is total lateness minutes (the count KPI is checked separately).

## Options for the paper (for the author to decide)

Pre-registered response: recentre on timing (RQ3) or the reversal condition (H2), or add heterogeneous intervention costs. In light of the diagnosis, the most promising options are:

1. **Forward-looking vs. backward-looking triage under capacity (recommended to explore first).** Core finding candidate: under customer time windows, the common "behind schedule" exception rule recovers only 4–44% of what predicted-risk triage recovers, because promise exposure is positional. Contribution: decomposing the value of prediction into *where* lateness will occur versus *how late* the route is, plus capacity–value curves showing that a 5% intervention budget suffices. Risk: may be read as "prediction beats no prediction". It needs the analytical piece (when delay-based and risk-based rankings coincide) to be publishable in TR-E.
2. **How much to change, not only which route (disruption-limited resequencing).** The oracle moves 25–64 stops, which is unrealistic. The decision becomes a joint choice of routes and resequencing depth (a multiple-choice knapsack on value per moved stop), where diminishing returns in depth could make "many light interventions" beat "few deep ones". Needs a new Gate (value curve Δ(k) for k = 1, 3, 5 relocations).
3. **Do-no-harm intervention under prediction noise.** Resequencing worsens 5–15% of routes even with a perfect model. A rule for *whether* to intervene on a route (an uncertainty-adjusted threshold) is a clean H2-type contribution, closest to the original Proposition 2.
4. **Drop the "recoverable vs. latest" headline.** Keep recoverable lateness only as a diagnostic; Proposition 1 becomes an explanation of why risk ranking works on real routes (comonotonicity at the top of the distribution).
5. **Keep the recoverability story only as a boundary condition.** The exploratory runs hint that recoverability ranking gains up to about 14% when lateness is widespread. A paper built on this would need a regime map (lateness concentration × capacity) and would rest on a modest effect; on current evidence it is the weakest option.

## Change log

| Date | Change |
|---|---|
| 2026-11-25 | Initial report |
