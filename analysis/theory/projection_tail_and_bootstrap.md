# Check 3c (tail form of Proposition 5) and bootstrap intervals: pre-specification

**Written 2026-11-27, before computing any quantity below.**

Like check 3b, this specification follows the gates but is written before its statistics are computed. It is therefore a weaker standard than a frozen gate.

**Known when this was written:**
- the cell-level efficiencies of L0–L4 and the captured-risk losses of L1 and L2 (Gates 2c and 3);
- the maximum form of the projection bound $B_{ii}=2K\max_r\tau_r/V(S_\rho)$, which is vacuous in every cell and has Spearman 0.56 with the L2 loss over the 28 W cells;
- that `tau_ii_min` is stored per route in `analysis/gate3/results/gate3_A.csv` and `gate3_B.csv`.

**Not known:** any value of $T_{2m_K}(\tau)$, and any bootstrap interval.

Code: `code/resequencing/bootstrap_tail_check.py` (written after this file). Outputs: `analysis/theory/tail_bootstrap/`.

---

## Part 1. Tail form of Proposition 5 (check 3c)

**Statistic.** For each cell:
- set $K=\max(1,\mathrm{round}(0.05N))$; since $K\le N/2$, $m_K=K$;
- with $\tau_r$ = `tau_ii_min` and $V(S_\rho)$ the sum of the $K$ largest `rho_min`, compute
$$B_{\rm tail}=T_{2K}(\tau)/V(S_\rho).$$

**Loss.** `loss_rho_L2` $=1-V(S_\pi)/V(S_\rho)$, exactly as in `analyze_gate3.mechanism`. Ties are broken by a stable sort in route order.

**Cell sets.**

| Set | Cells | Role |
|---|---|---|
| P (primary) | the 28 W cells, Gate 3 arm A at $\sigma_U=0$ | the set on which the maximum form gave Spearman 0.56 |
| S1 | the 40 arm-B cells | secondary |
| S2 | the 36 arm-A cells with $\sigma_U>0$ | exploratory |

**Must hold in every cell, because they are theorems:**
- $B_{\rm tail}\ge$ `loss_rho_L2`;
- $B_{\rm tail}\le B_{ii}$.

If either fails, the code is wrong: stop and debug.

**What is reported:**
- the share of cells with $B_{\rm tail}<1$, overall and split by $\omega\ge0.5$ and $\omega<0.5$;
- Spearman($B_{\rm tail}$, `loss_rho_L2`) in P and in S1;
- the median ratio $B_{\rm tail}$/loss over cells with loss $>0.01$.

**Expectation.** Remark 1 gives the direction:
- $B_{\rm tail}$ is below 1 in most cells with $\omega\ge0.5$;
- it is at or above 1 in most ETA and DAY cells;
- Spearman is probably above the 0.56 of the maximum form. Whether it reaches 0.80 is open, because the projection's loss is largest exactly where $\tau$ is large.

**Wording rule.**

| Result | Wording |
|---|---|
| Spearman $\ge0.80$ in both P and S1, and $B_{\rm tail}<1$ in at least half of P | Section 5.3 says the tail form is informative and ranks cells by the projection's loss, with the numbers |
| Spearman $\ge0.80$ in both, but $B_{\rm tail}\ge1$ in more than half of P | "ranks the cells but is mostly vacuous" |
| Otherwise | Report the numbers in one sentence; Proposition 5 stays a guarantee |

In every case, the sentence "has not been evaluated" is replaced by the result.

---

## Part 2. Bootstrap intervals for simulated efficiencies

**Purpose.** The intervals show the sampling uncertainty that comes from which routes were simulated. Each route's simulated $\rho_r$, $\Delta_r$ and scores are treated as fixed, so Monte Carlo error within a route is not included.

**Method.**
- **Resamples.** $B=2{,}000$, with percentile 95% intervals and seed 2032.
- **Efficiency.** Efficiency is the ratio of sums in equation (10). For each resample:
  - set $K=\max(1,\mathrm{round}(0.05N^*))$;
  - recompute every selection by a stable sort of the scores;
  - sort the resampled rows by original route order, so that ties are broken as in the point estimate.
- **Gate 2c.** Resample routes (2,999). Use one resample of route indices for all 28 cells, so that cells are compared on the same routes.
- **Gate 3.** Resample whole station-days (170) as a cluster bootstrap, because routes of a station-day share a date and a station.

**Quantities.**

| Gate | Quantity | Cells |
|---|---|---|
| 2c | Eff(L1), Eff(L2), Eff(L3), and Eff(L2) − Eff(L1) | all 28 cells, $K/N=5\%$ |
| 2c | Eff(L1) against 0.70 (primary criterion C1) | MIXB, 35% check-in, base noise |
| 3 | Eff(L2) − Eff(L1) (H5a) | the five cells at $\sigma_U=0.5$ with $\omega\ge0.5$ |
| 3 | Eff(L3) against 0.80 (H5b) | the 12 cells at $\sigma_U=0.5$ |
| 3 | Eff(L1) against 0.80 (H3-D) | arm B, MIXB, base, $\varphi=0.75$ |

**Expectation.** With about 150 selected routes per cell and lateness concentrated in few routes:
- half-widths of a few hundredths, roughly 0.02–0.06;
- the C1 interval excludes 0.70;
- the H5a intervals exclude 0;
- the intervals of the two H5b cells that fell below 0.80 (0.786 and 0.792) include 0.80.

**Wording rule.**

| Result | Wording |
|---|---|
| C1: upper end below 0.70 | "outside sampling noise" |
| C1: otherwise | The claim is qualified |
| H5a: all five intervals exclude 0 | Stated as such |
| H5a: otherwise | The cells whose interval includes 0 are named |
| H5b: a failing cell's interval includes 0.80 | The text says it fell below 0.80 "by less than the sampling uncertainty". The pre-specified outcome (H5b failed) is unchanged. |
| Table 3 | A note gives the largest 95% half-width across cells for each level, and the full intervals go to the replication files |

If a difference between levels that the text claims has an interval that includes 0, the text names it.

No gate decision is revisited.

---

## Results (computed 2026-11-27, after this specification was written)

Outputs are in `analysis/theory/tail_bootstrap/`: `check3c_cells.csv`, `bootstrap_gate2c.csv`, `bootstrap_gate3.csv`, `summary.json`, `run.log`.

### Part 1 (check 3c)

Both theorem checks hold in every cell: $B_{\rm tail}\ge$ loss, and $B_{\rm tail}\le B_{ii}$.

| Set | Spearman (tail) | Spearman (maximum form) | Share below 1 | … if $\omega\ge0.5$ | … if $\omega<0.5$ | Median bound/loss |
|---|---|---|---|---|---|---|
| P (28) | 0.877 | 0.557 | 0.36 | 0.80 | 0.11 | 8.3 |
| S1 (40) | 0.849 | 0.901 | 0.35 | 0.64 | 0.00 | 9.9 |
| S2 (36, exploratory) | 0.810 | 0.662 | 0.53 | 0.71 | 0.41 | 8.0 |

**Wording rule applied:** Spearman is at least 0.80 in both P and S1, but $B_{\rm tail}\ge1$ in more than half of P. The text therefore says the tail form "ranks the cells but is mostly vacuous", and Proposition 5 stays a guarantee.

**Against the expectation:**
- direction as expected: below 1 mostly where $\omega\ge0.5$, and rarely otherwise;
- Spearman above 0.56 as expected, and it reached 0.80.

### Part 2 (bootstrap)

| Criterion | Value | 95% interval | Reading |
|---|---|---|---|
| C1, Gate 2c | 0.603 | [0.501, 0.689] | outside sampling noise |
| H5a, five cells | L2 − L1 | every interval excludes 0 | |
| H5b, MIXB | 0.786 | [0.720, 0.851] | includes 0.80 |
| H5b, DAY | 0.792 | [0.731, 0.854] | includes 0.80 |
| H3-D | 0.895 | [0.838, 0.931] | |

The H5b wording is "smaller than the sampling uncertainty". The pre-specified outcome is unchanged.

**Largest Gate 2c half-widths:**

| Quantity | Half-width |
|---|---|
| L1 | 0.099 |
| L2 | 0.183 |
| L3 | 0.046 |
| L2 − L1 | 0.169 |

**Against the expectation.** These half-widths are wider than expected (0.02–0.06). The widest L2 intervals are in the ETA and DAY cells, where few routes have positive projected lateness.

Every L1 versus L2 difference that the text discusses has an interval that excludes zero:
- ETA and DAY at base noise: L2 < L1;
- MIXB, MIXB20, MIXB50 and TW at base noise: L2 > L1.
