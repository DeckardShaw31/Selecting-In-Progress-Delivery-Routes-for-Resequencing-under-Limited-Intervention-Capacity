# Gate 2 pre-registration: behind schedule vs. at risk (option 1)

Written 2026-11-25, after Gate 1 and **before** any Gate 2 simulation. Gate 1 results that motivate this gate are known (`analysis/gate1/gate1_report.md`); in particular, on the Gate 1 cells the "delay so far" rule captured 78–99% (ETA promises) and 4–44% (customer windows) of the oracle value. Changes after seeing Gate 2 results must be logged at the bottom.

## Re-framed question

A dispatcher who can intervene on only K of N routes usually triages by **how far behind schedule** a route is (deviation-based exception management). The alternative is to triage by **how much promised-time lateness the route is about to incur** (forward-looking risk). Gate 2 asks:

1. **Gap:** under realistic promise structures, how much of the lateness that forward-looking triage removes does behind-schedule triage remove?
2. **Decomposition (information ladder):** how much of the gap closes with (L2) a deterministic projection of the current delay onto the remaining promises (no forecasting model), (L3) adding generic travel/service-time variability (no route-specific learning), and (L4) adding route-specific learning from progress so far?
3. **Mechanism (positional exposure):** is the efficiency of behind-schedule triage governed by how unevenly lateness exposure is spread across routes?

## Mechanism and stylized prediction

Write a route's expected lateness as a function of its current delay D: rho_r = E_r(D_r), where the exposure function E_r depends on where the route's remaining promises sit (their planned slack). Its slope b_r = dE_r/dD equals the expected number of late promises ("exposure mass"; derivative of E[(D + xi - s)^+] is P(D + xi > s)). If all routes shared the same E_r, ranking by D would equal ranking by risk. In the linearized model rho_r = b_r * D_r with b_r independent of D_r and E[D] = 0, corr(rho, D) = 1 / sqrt(1 + CV(b)^2), so the efficiency of behind-schedule triage should fall as the coefficient of variation of exposure mass across routes grows. This is a stylized prediction to be tested, not an assumption of the experiment.

## Data and execution model

Identical to Gate 1 (Amazon evaluation split, all 3,052 routes; same lognormal model, posterior update, SAA resequencing oracle with 48 in-sample and 256 out-of-sample scenarios, seeds and common random numbers). Noise levels: **base** and **high** (Gate 1 values). Check-in at f ∈ {0.35, 0.50}. Route-level service shock off (as pre-registered for Gate 1).

## Promise regimes

- **ETA:** every stop promised by planned ETA + 60 min (Gate 1).
- **TW:** Amazon customer windows only; stops without windows carry no promise (Gate 1).
- **MIXB (primary realistic regime): day-definite + windows.** Stops with a customer window are promised by that window (waiting if early); all other stops are promised by the planned route end + 60 min (one day-definite deadline per route). All stops weight 1. Represents standard parcel practice: most parcels due by end of day, some time-windowed.
- **MIXA (secondary): communicated ETA + windows.** Every stop promised by ETA + 60 min; windowed stops by the earlier of their window end and ETA + 60 min (waiting if early).

## Information ladder (selection scores, all computed at the check-in)

| Level | Score | Information used |
|---|---|---|
| L0 | random | none |
| L1 | delay so far D_r (realized completion time at check-in minus planned) | backward-looking deviation only |
| L2 | lateness of the **deterministic projection**: remaining stops in current order, mean travel and service times, starting from the realized position and time | current delay + positions of remaining promises; no forecasting model |
| L3 | expected lateness under the **prior** travel/service model (route shock not updated) | + generic variability |
| L4 | expected lateness under the **posterior** (route shock learned from completed legs) = Gate 1 risk | + route-specific learning |
| ORC | true recoverable lateness Delta (Gate 1 oracle) | upper bound |

Value of a selection = sum of Delta (out-of-sample expected lateness removed by resequencing) over the selected routes, as in Gate 1. Efficiency of a level = its value / value of L4 at the same K.

## Metrics

- V_L1/V_L4, V_L2/V_L4, V_L3/V_L4 at K/N ∈ {2%, 5%, 10%}.
- Share of the L1→L4 gap closed by L2 and by L3: (V_L2 − V_L1)/(V_L4 − V_L1), (V_L3 − V_L1)/(V_L4 − V_L1).
- Exposure mass b_r = expected number of late promises under the posterior (current sequence); CV(b) across routes with b_r > 0 and across all routes; corr(D, rho) (Pearson, across routes).

## Decision rule (thresholds set now)

- **G2a (practical gap, primary regime):** at f = 0.35, base noise, K/N = 5%: V_L1/V_L4 ≤ 0.70 in **MIXB**.
- **G2b (mechanism):** across the 16 cells (4 regimes × 2 check-ins × 2 noise levels) at K/N = 5%: (i) Spearman(V_L1/V_L4, corr(D, rho)) ≥ 0.70 and (ii) Spearman(corr(D, rho), 1/sqrt(1 + CV(b)^2)) ≥ 0.70, with CV(b) computed over all routes.

| Outcome | Action |
|---|---|
| G2a and G2b pass | **GO** with option 1: "behind schedule is not at risk" with positional exposure as the mechanism. |
| G2a fails, G2b passes | **NARROW**: the gap is real only where promises are heterogeneous (windowed deliveries); scope the paper to attended/windowed delivery or present a regime map. |
| G2b fails | **REVISE**: the gap exists but the exposure mechanism does not explain it; the paper would be empirical only, which is weaker for TR-E. |

The decomposition (L2, L3 shares) is reported without a threshold. Pre-stated interpretation: if L2 closes ≥ 80% of the gap, the managerial message is "promise-aware projection suffices, no forecasting model needed"; if ≤ 50%, "forecasting matters"; in between, both.

## Change log

| Date | Change | Reason |
|---|---|---|
| 2026-11-25 | Initial pre-registration | — |
| 2026-11-25 | Bug fix before the main run: MIXA ETAs now use the planned schedule including waiting for windows to open (previously the no-wait ETA, which promised some windowed stops before their window opened). Delay so far (L1) is measured against the regime-consistent planned schedule (with waiting where the regime waits); Gate 1's TW "delay so far" used the no-wait plan. | Found in a 32-route smoke test (implementation check only; no Gate 2 results inspected beyond regime means). |
| 2026-11-25 | After the main run: exploratory analysis of the failed mechanism test (corr(b, D), rank correlations) added in `exploratory_mechanism.csv`. Not part of the decision; any re-specified mechanism must be confirmed on independent data. | Decision rule applied first: outcome REVISE. |
