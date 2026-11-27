# Gate 1 pre-registration: is recoverable lateness material and different from predicted lateness?

Written 2026-11-25, **before** any route was simulated. Changes made after seeing results must be logged at the bottom with a reason.

## Question

On real last-mile routes, if a dispatcher resequences the remaining stops of a route at a mid-route check-in:

1. **Materiality:** what share of the route's expected lateness can resequencing remove (recoverable share)?
2. **Divergence:** does ranking routes by expected lateness (risk) select different routes from ranking by recoverable lateness, and how much does that cost at small capacity?

Gate 1 is an **upper-bound check** under a perfectly specified predictive model (the dispatcher knows the data-generating process and updates the route shock from observed progress). It does **not** test H1 under realistic prediction error; that is Experiment 2.

## Data

- 2021 Amazon Last Mile Routing Research Challenge, **evaluation split** (3,072 routes; `data/raw/amazon_lmrrc_eval/`).
- Observed (O): stop sets, actual stop sequences, departure date and time, package time windows, planned package service times, average realized travel times between stops.
- Derived (D): stop service time = sum of planned service times of packages at the stop; stop deadline = earliest package time-window end at the stop; stop earliest start = latest package time-window start at the stop; actual sequence treated as the dispatched plan.
- Assumed (A): all stochastic parameters below.

## Execution model (A)

- Travel time on leg (i, j): `T_ij = tbar_ij * exp(theta_r + eps_ij - (sR^2 + sE^2)/2)` with route-day shock `theta_r ~ N(0, sR^2)` and leg noise `eps_ij ~ N(0, sE^2)`. Means are preserved.
- Service time: `S_i = sbar_i * exp(zeta_i - sS^2/2)`, `zeta_i ~ N(0, sS^2)`.
- Noise levels (sR, sE, sS): **low** (0.10, 0.15, 0.25), **base** (0.15, 0.25, 0.40), **high** (0.25, 0.35, 0.60).
- Information at check-in: the realized times of completed legs and services. The route shock is updated by the normal–normal posterior using the observed log-ratios of completed legs (`log(T/tbar) + (sR^2+sE^2)/2 = theta + eps`). Future legs and services are drawn from the posterior predictive.

## Promise regimes

- **TW (observed windows):** lateness at stops with an Amazon time window, `(A_i - deadline_i)^+`. Early arrival waits until the window opens. Stops without a window carry weight 0 but must still be visited.
- **ETA (communicated windows, derived):** at dispatch, planned arrival `ETA_i` is computed along the actual sequence with mean times; promise = `ETA_i + h`, with h = 60 minutes (a 2-hour window centred on the ETA). Every stop carries weight 1. No waiting for early arrival. Assumed to mimic carriers that send ETA windows on the morning of delivery.

## Check-in and intervention

- Check-in after a fraction f of stops has been served: f ∈ {0.35, 0.50}.
- Intervention = any reordering of the remaining stops (the return to the station is not penalized).
- Resequencing oracle: SAA local search (Or-opt with segment lengths 1–3, plus 2-opt) on S = 48 in-sample scenarios, started from both the current sequence and the earliest-deadline-first sequence; keep the better solution.
- **All reported values are out of sample**: rho (current sequence), lambda (resequenced), Delta = rho − lambda, each on 256 fresh scenarios from the same posterior predictive (common random numbers for both sequences).

## Metrics

- Recoverable share `R = sum(Delta) / sum(rho)` over routes with rho > 0.
- Distribution of Delta and of `Delta/rho` across routes (median, P90, share of routes with Delta/rho ≥ 0.25).
- Spearman correlation between rho and Delta.
- For K/N ∈ {5%, 10%, 20%} of routes in the sample: value captured by risk top-K (`sum Delta` over the top-K by rho) as a fraction of the value of recoverability top-K (`V_risk / V_rec`), and the overlap of the two top-K sets.

## Decision rule (thresholds set before running)

Evaluated at f = 0.35, base noise:

- **G1a (materiality):** R ≥ 0.10 in at least one promise regime.
- **G1b (divergence):** at K/N = 10%, `V_risk / V_rec ≤ 0.85` in the same regime.

| Outcome | Action |
|---|---|
| G1a and G1b hold | **Go**: proceed with the paper as framed. |
| G1a holds, G1b fails | **Revise**: risk ranking is near-optimal on real routes; recentre the paper on timing (RQ3) and the reversal condition (H2), or add heterogeneous intervention costs. |
| G1a fails | **Revise (add a lever)**: intra-route resequencing removes too little; add stop transfer to a neighbouring route or a rescue vehicle before going further. |

These thresholds are the analyst's proposal; the author may override them, but any override must be logged below.

## Sample

Routes from the evaluation split with at least 10 stops. Base configuration on all eligible routes if runtime allows, otherwise a random sample of at least 600 routes (seed 2026) stratified by station; the low and high noise levels on the same sample.

## Change log

| Date | Change | Reason |
|---|---|---|
| 2026-11-25 | Initial pre-registration | — |
| 2026-11-25 | Data note: the evaluation split contains 3,052 routes (the registry page said 3,072). All 3,052 have at least 10 stops, so the main run uses all of them instead of a sample. | Observed after download, before any simulation. |
| 2026-11-25 | **Exploratory additions (not part of the decision rule):** (a) ETA promise half-width 30 min (1-hour windows); (b) a route-level service-time shock `phi_r ~ N(0, sP^2)`, sP = 0.20, updated from observed service times like the travel shock. | Added after a 32-route timing pilot showed that expected lateness is rare under the pre-registered settings (about 9% of routes late). The Gate 1 decision uses only the pre-registered configuration; the additions are reported separately and labelled exploratory. |
