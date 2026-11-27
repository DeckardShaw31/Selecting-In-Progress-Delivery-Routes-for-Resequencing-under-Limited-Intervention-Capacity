# Gate 3 pre-registration: imperfect prediction, delay structure, within-day selection and the loss bound

Written 2026-11-26, **before any Gate 3 simulation was run**. Three things were done before writing, all descriptive (no triage score and no lateness outcome was computed):

1. Calibration facts were measured (`code/resequencing/gate3_calibrate.py` → `analysis/gate3/calibration.json`).
2. The route sample was drawn and frozen (`code/resequencing/gate3_sample.py` → `analysis/gate3/sample_routes.csv`).
3. The Proposition 2 terms have **not** been computed on any earlier gate's results.

Any change made after seeing results must be logged at the bottom with a reason.

## Purpose

Gates 2 and 2c established the regime map under three simplifications. The predictor was perfectly specified, the noise levels were assumed, and selection was pooled across station-days. Framework Section 5.4 lists these as the ways the claim can still be wrong. Gate 3 tests each of them, and also tests the analytical claim:

| Question | Experiment | Hypothesis |
|---|---|---|
| Does promise-aware triage keep its edge when the carrier's predictions are wrong? | E6 (Arm A) | H5 (primary: decides GO / GO-NARROW / REVISE) |
| Does a larger route-wide share of delay, at equal total variance, shrink the behind-schedule loss? How large is the loss at the calibrated delay structure? | E7 (Arm B) | H3, H3-D |
| Does the Proposition 2(i) term rank cells by behind-schedule loss? | E5 | H4 |
| Does the gap survive when each station-day selects from its own routes? | Part of E8, moved here | H6 |

## Calibration facts (measured 2026-11-26, before this document)

**Station-days in the Amazon data.**
- Training split: 6,112 routes in 521 station-days (17 stations, 39 dates). The median station-day has 10 routes (IQR 5–15, maximum 56), and 80% of routes sit in station-days with at least 10 routes.
- Evaluation split: median 2 routes per station-day (maximum 25).
- So within-day selection can only be studied on the training split.

**Delay structure in observed data (Konovalenko et al. calibration weeks).**
- Data: 5,304 calibration routes; 75,420 legs from 5,290 of them. Legs with 0 < inter-arrival time < 99th percentile (152.9 min) and a non-depot destination.
- The residual log(actual inter-arrival / fleet-pace prediction) was split into a route effect and a leg residual, using the fleet pace of Gate 2d (19.3 min + 1.67 min/km).
  - Route effect: τ = 0.49.
  - Leg residual: κ = 1.32.
  - Route share τ²/(τ² + κ²) = **0.12**, or 0.10 on the 5,014 routes whose arrival times are monotone.

**The same statistic in the simulator.** Computed on 300 training routes, measured against the carrier's mean times.
- At the assumed noise levels, the route share is 0.04 (low), 0.04 (base) and 0.05 (high).
- At equal total travel log-variance, put a share φ of that variance on the route-wide shock:
  - φ = 0.50 gives a route share of 0.07–0.08;
  - **φ = 0.75 gives 0.109 (base total) and 0.116 (high total), matching the observed share**;
  - φ = 0.90 gives 0.13–0.14.
- The assumed levels correspond to φ = 0.265 (base) and 0.338 (high).

**Reading.**
- The assumed noise levels probably understate how much of the delay is route-wide. So the regime map is re-measured at φ = 0.75, labelled **D (calibrated)**.
- Both observed components are measured against a crude predictor, so the observed share is approximate. Route-level predictor error inflates τ, and leg-level error inflates κ.
- The observed leg residual (κ = 1.32, against 0.29 in the simulator at base noise) reflects how crude the fleet-pace predictor is: every leg gets the same intercept and there is no stop-level service time. It is therefore **not** used to set noise levels or prediction error.

## Data and sample

- **Source:** Amazon LMRRC training split, preprocessed as in Gate 2c (`amazon_io.py`, `data/processed/train_routes.pkl`).
- **Sample: whole station-days.**
  - Station-days with at least 10 routes (each route with at least 10 stops) were shuffled with seed 2029 and taken in order until the sample held at least 3,000 routes.
  - Result: **3,028 routes in 170 station-days** (15 stations, 38 dates; 10–55 routes per station-day, median 15; 148.6 stops per route on average; 6.8% of stops with a customer window).
  - 1,469 of these routes were also in the Gate 2c sample. The overlap is acceptable because Gate 3 tests new hypotheses, and Arm A's perfectly specified setting (W) gives a paired baseline on the same routes.
- **Why whole station-days:** within-day selection (H6) can then be evaluated on the same runs, at no extra simulation cost.

## Common design (as in Gate 2c unless stated)

- **Execution model.** Lognormal travel times with a route-wide shock θ and leg noise ε, mean-preserving. Service time is the planned time × exp(ζ − σ_S²/2). There is no route-wide service shock.
- **Check-in:** after a share f of stops.
- **Promise regimes:** ETA, MIXA, DAY, MIXB, MIXB20, MIXB50 and TW, defined as in Gate 2c. The ETA half-width is 60 min. The day-definite deadline is the planned end (with waiting) + 60 min. Synthetic windows are drawn from the frozen evaluation-split offsets.
- **Resequencing oracle:** SAA with 48 in-sample and 256 out-of-sample scenarios, Or-opt + 2-opt, objective in minutes.
- **Randomness:** common random numbers across scores and settings. Simulation seed 2029.
- **Value and efficiency.**
  - Value: V(S) = sum of out-of-sample recoverable lateness Δ over the top-K routes ranked by score S (ties broken by route order).
  - Efficiency: V(S) / V(L4*), where L4* is risk triage under the true law.
  - Selection is pooled across the sample at K/N = 5% (primary); 2% and 10% are also reported.
- **Ladder.**
  - L0: random.
  - L1: delay so far against the regime-consistent plan (`delay_plan_min`).
  - L2: deterministic projection of the remaining promises with mean times, from the current time and position.
  - L3: prior forecast (route shock not learned).
  - L4: posterior forecast.
  - Oracle: recoverable lateness itself.
- **ω (share of a cell's expected late promises that are customer windows):** ω = Σ expected late window promises / Σ expected late promises, out of sample under the true law. As in Gate 2c, under ETA promises windows are not promises, so ω = 0 there.

## Arm A: imperfect prediction (E6; H5)

**Setting.**
- **The carrier's model.** The Amazon mean travel times (tbar) and planned service times (sbar) are the carrier's predictions. They define:
  - the plan;
  - the ETA and day-definite promises;
  - the anchors of the synthetic windows;
  - L1;
  - every forecast.
- **The truth.** True mean times differ from the carrier's by persistent multiplicative errors:
  - travel on arc (i, j): tbar_ij · exp(σ_U z_ij − σ_U²/2);
  - service at stop i: sbar_i · exp(σ_U z′_i − σ_U²/2).
- **The error draws.** The z are independent standard normals, drawn once per route from a route-specific stream. The same z are used for every σ_U, regime, check-in and noise level.
- **Properties of the errors.** They are unbiased in level. They are independent across arcs, which is the hardest case for learning: completed legs say nothing about future arcs except through the route-shock posterior.
- **Error levels σ_U:**
  - 0, labelled W (perfectly specified; the Gate 2c design on the new sample);
  - 0.25 (good): a mean absolute error of about 20% per leg mean;
  - **0.50 (poor; primary)**: about 39%;
  - 1.00 (stress): about 77%.
- **Reference:** the simulated leg residual SD against the carrier's means at base noise is 0.29 at σ_U = 0, 0.35 at 0.25 and 0.49 at 0.50.

**What each party sees.**
- **Truth and oracle.**
  - Execution uses the true means and the same random draws as W.
  - The oracle and L4* use the true law. The posterior of θ is computed from completed legs relative to the true means.
  - Δ, ω and the late-promise counts are computed out of sample under the true law.
- **Dispatcher (forecasts).** All use the carrier's means:
  - L2: deterministic projection from the true current time and position;
  - L3: prior forecast;
  - L4f: posterior forecast, with the θ posterior computed from completed legs relative to the carrier's means, so it absorbs part of the error.
- **Additional variant at σ_U = 0:** the dispatcher's dispersion parameters (σ_R, σ_E, σ_S) are all scaled by 0.5 or by 2. This gives two extra versions of L3 and L4f, reported without thresholds.

**Cells.**
- W: 7 regimes × f ∈ {0.35, 0.50} × noise ∈ {base, high} = 28 cells.
- Each σ_U > 0: regimes ETA, DAY, MIXB, MIXB20, MIXB50 and TW × f = 0.35 × noise ∈ {base, high} = 12 cells.

## Arm B: delay structure at equal total variance (E7; H3, H3-D)

- **Model:** perfectly specified (σ_U = 0).
- **Travel noise:** (σ_R, σ_E) = (√(φV), √((1 − φ)V)), where V = σ_R² + σ_E² of the base level (0.085) or the high level (0.185). Service σ_S keeps that level's value (0.40 or 0.60). The posterior uses the cell's true (σ_R, σ_E).
- **Route-wide share φ:** {0.10, 0.25, 0.50, **0.75 (D)**, 0.90}.
- **Regimes:** DAY, MIXB, MIXB50, TW. **Check-in:** f = 0.35.
- **Cells:** 2 × 5 × 4 = 40.

## Mechanism statistics (E5; H4)

**Inputs, per cell:**
- risk ρ_r = L4*, the true posterior expected weighted lateness of the current sequence, out of sample;
- delay D_r = the L1 score;
- K = round(0.05 N).

**Primary statistic: the Proposition 2(i) term.**
- ε* = min over nondecreasing ḡ of max_r |ρ_r − ḡ(D_r)| = ½ · max over route pairs with D_a ≤ D_b of (ρ_a − ρ_b)⁺. This is the error of the L∞ isotonic fit of risk on delay.
- Relative bound: B_i = 2K ε* / V_ρ(S_ρ), where V_ρ(S) = Σ_{r∈S} ρ_r.
- Proposition 2(i) implies 1 − V_ρ(S_D) / V_ρ(S_ρ) ≤ B_i in every cell. This is checked in code; any violation is a bug and must be fixed before reporting.

**Secondary statistics (no thresholds):**
- (a) H_iso = 1 − R² of the least-squares isotonic regression of ρ on D (an average-type measure of heterogeneity).
- (b) ω.
- (c) The Proposition 2(ii) term for the projection.
  - Per route: τ_r = Σ_i w_ri · Ê[|A_ri − Â_ri| · 1{max(A_ri, Â_ri) > d_ri}], where A is the out-of-sample service-start time under the true law and Â is the projected one.
  - Per cell: B_ii = 2K · max_r τ_r / V_ρ(S_ρ).
  - B_ii is compared with the projection's loss in captured risk and checked to be a valid bound.
- (d) The share of cells in which the bound is vacuous (B ≥ 1).

**Replication:** the same statistics are computed on the Gate 2c route results (28 cells), after this document is frozen.

## Within-station-day selection (H6; part of E8 moved here)

- **Capacity:** station-day d with N_d routes gets K_d = max(1, round(0.10 N_d)).
- **Value:** V_day(S) = Σ_d of Δ over the top-K_d routes of station-day d ranked by S.
- **Efficiency:** V_day(S) / V_day(L4*).
- **Cells:** computed on the 28 W cells. K_d = 1 is also reported.

## Hypotheses and decision rule (thresholds set now)

**Primary: H5 decides the outcome.** Both parts are tested at σ_U = 0.50, f = 0.35 and K/N = 5%.

- **H5a (projection under imperfect prediction):** in every cell where ω ≥ 0.50, V(L2) ≥ V(L1).
  - Candidate cells: MIXB, MIXB20, MIXB50 and TW, at base and high noise.
  - ω is computed under the true law of each cell.
  - "Qualifying cells" below means the cells with ω ≥ 0.50.
- **H5b (prior forecast under imperfect prediction):** V(L3) / V(L4*) ≥ 0.80 in all 12 cells.

| H5a | H5b | Outcome |
|---|---|---|
| holds in all qualifying cells | pass | **GO**: the regime map and ladder claims stand under imperfect prediction; proceed to Phase 5. |
| holds in all qualifying cells | fail | **GO-NARROW (forecast)**: the projection claim stands. The practice rule "a prior forecast is within about 12% of risk triage everywhere" is restated as conditional on prediction quality. |
| fails in at most one-third of qualifying cells | either | **GO-NARROW (projection)**: report the cells where behind schedule wins under imperfect prediction; the practice rule gains a prediction-quality condition. |
| fails in more than one-third of qualifying cells | either | **REVISE**: promise-aware triage pays only with accurate predictions. Prediction quality becomes a third axis of the regime map, and the abstract and contributions are rewritten around it before Phase 5. |

**Secondary hypotheses.** These decide claim wording only.
- **H3 (delay structure):** tested in four series (MIXB50 and TW, at base and at high total variance), with K/N = 5%.
  - Condition 1: L1 efficiency at φ = 0.90 minus L1 efficiency at φ = 0.10 ≥ 0.10.
  - Condition 2: Spearman(φ, L1 efficiency) over the five φ levels ≥ 0.90, which allows at most one adjacent swap.
  - H3 is supported if at least 3 of the 4 series meet both conditions.
- **H3-D (gap at the observed window share under the calibrated delay structure):** at φ = 0.75, base total variance, MIXB, f = 0.35 and K/N = 5%, L1 efficiency ≤ 0.80.
- **H4 (loss bound as a mechanism statistic):** across the 28 W cells, Spearman(B_i, captured-risk loss of L1) ≥ 0.80. The captured-risk loss is 1 − V_ρ(S_D)/V_ρ(S_ρ). The same correlation with the recoverable-lateness loss (1 − L1 efficiency) is also reported.
- **H6 (within station-days):** in MIXB at f = 0.35 and base noise, within-day L1 efficiency ≤ 0.80.

| Result | Consequence for the paper |
|---|---|
| H3 supported | Delay structure is stated as the second axis of the regime map. A heat map (window share × φ) is a main exhibit. |
| H3 not supported | The attenuation under high noise is reported as observed, without a causal interpretation in terms of the route-wide share. |
| H3-D met | The paper may state that behind-schedule triage forgoes at least 20% at observed window shares under the calibrated delay structure. |
| H3-D not met | The loss at observed window shares is given as a range from the φ = 0.75 value to the base value. The abstract leads with window-dominated settings, not with the observed share. |
| H4 met | The Proposition 2(i) term is reported as the mechanism statistic alongside ω. |
| H4 not met | The bound is presented as a guarantee only; ω remains the empirical mechanism statistic. |
| H6 met | The gap is stated to persist within station-days. |
| H6 not met | Within-day numbers lead, and pooled numbers are labelled as an upper bound on the gap. |

**Reported without thresholds:**
- all ladder efficiencies per cell at 2%, 5% and 10%, with random selection and L4*/oracle;
- the cost of misspecification, V(L4f)/V(L4*);
- the σ_U curve (0, 0.25, 0.50, 1.00);
- the variants with misspecified dispersion;
- the W replication of Gate 2c (Gate 2c: MIXB, f = 0.35, base, L1 efficiency 0.603);
- the Arm B heat map for every regime and ladder score;
- B_ii and the other secondary statistics;
- within-day results for all W cells;
- run time.

## Implementation and verification (written after this document, before the main run)

- **Scripts.**
  - `code/resequencing/gate3.py` reuses the Gate 1–2c helpers and adds:
    - the true-mean perturbation;
    - the separation between truth and dispatcher;
    - the split of travel noise into route-wide and leg shares (Arm B);
    - a per-stop service-start function for B_ii.
  - `code/resequencing/analyze_gate3.py` applies the decision rules above.
  - Both are written before the main run. The decision is recorded before any exploratory analysis.
- **Verification before the main run:**
  1. With σ_U = 0, the base split and seed 2027, `gate3.py` reproduces the `gate1.py` Gate 2c row values on 20 routes.
  2. Lateness computed from the per-stop service-start times equals `eval_seq` on random instances.
  3. The closed form for ε* matches a brute-force L∞ isotonic fit on random data.
  4. Both Proposition 2 bounds hold in every cell (checked in the analysis).
- **Compute:** about 104 route-configurations per route × 3,028 routes. At Gate 2c's rate, about 3–4 hours on 14 cores.

## Not covered by Gate 3

- A predictor learned on observed data. Gate 2d's fleet-pace projection is the observed counterpart (κ = 1.32, far cruder than σ_U = 1.00), and E9 extends it.
- A dispatcher who also resequences with the wrong model. Value is always measured with oracle resequencing, so that triage information is isolated.
- Move caps, adoption and intervention cost. These remain E8.

## Change log

| Date | Change | Reason |
|---|---|---|
| 2026-11-26 | Initial pre-registration | — |
