# Selecting In-Progress Delivery Routes for Resequencing under Limited Intervention Capacity: replication package

**Paper.** Nguyễn Nhật Thanh (International University, Vietnam National University, Ho Chi Minh City, Vietnam), "Selecting In-Progress Delivery Routes for Resequencing under Limited Intervention Capacity", submitted to the *European Journal of Operational Research*.

**Repository.** <https://github.com/DeckardShaw31/Selecting-In-Progress-Delivery-Routes-for-Resequencing-under-Limited-Intervention-Capacity>

**Archive.** [TO FILL after the first release: Zenodo DOI]

This package holds:
- the code: simulator, resequencing oracle, gate and analysis scripts, and numerical checks of the theory;
- the pre-specification documents of every experiment ("gate") and their change logs;
- the derived results behind every number, table and figure in the paper and its Online Supplement.

The raw data are **not** included. `data/README.md` says where to get them.

## Contents

| Path | What it is |
|---|---|
| `code/resequencing/` | Python code: simulator (`sim.py`, `gate1.py`), data readers, one script per gate, analysis scripts, and the checks of the theory (`theory_checks.py`, `test_sim.py`) |
| `analysis/gate1/` … `analysis/gate5/` | For each gate: its `preregistration.md` (pre-specification and change log), its report with the exact commands and seeds, and its route- and cell-level results and figures |
| `analysis/theory/` | Proof audit (`proof_audit.md`, which uses an older result numbering, mapped below); the specifications and outputs of check 3b (`average_bound_check.md`, `avg_bound/`) and check 3c with the bootstrap intervals (`projection_tail_and_bootstrap.md`, `tail_bootstrap/`); the route-shock checks (`route_shock/`) |
| `figures/` | The four figures of the paper (vector PDF, and 600 dpi PNG) |
| `data/README.md` | Where to put the raw data, and how to build the processed files |

The largest files are the route-level results in `analysis/*/results/` (up to 82 MB each).

## Software

- **Python 3.13.11**, with the package versions in `requirements.txt` (`pip install -r requirements.txt`).
- **numba** compiles the simulator on first use.
- Run every script from the repository root as `python -I code/resequencing/<script>.py …`; each script states its usage in its docstring.
- Each gate's report (`analysis/<gate>/*_report.md`) records the exact commands, seeds and run times.
- The full simulation gates take hours on a desktop computer. Their results are stored at route level, so the analysis scripts (`analyze_gate*.py`, `avg_bound.py`, `bootstrap_tail_check.py`, `gate5.py`) recompute the simulation tables and figures in minutes, without re-simulation.
- The observed-data analyses (Gates 2d and 4) start from the raw data.

## From the paper to the code

Section, table and figure numbers follow the journal version. "S" numbers refer to its Online Supplement.

| Paper item | Gate or check | Scripts |
|---|---|---|
| Triage target (§2, Remark 3) | Gate 1 | `gate1.py`, `analyze_gate1.py` |
| Primary criterion on the evaluation split (§5.1) | Gate 2 | `gate1.py`, `analyze_gate2.py` |
| Table 3, Figure 2, §5.1–5.3 | Gate 2c (confirmatory) | `gate1.py`, `analyze_gate2c.py --figures` |
| Table 4 (delivery routes), §6 | Gate 2d | `konovalenko_io.py`, `gate2d.py`, `gate2d_robustness.py` |
| §5.3–5.4, Figures 3–4 | Gate 3 (robustness) | `gate3_calibrate.py`, `gate3_sample.py`, `gate3_verify.py`, `gate3.py`, `analyze_gate3.py --figures` |
| Tail-average bound (§5.3) | Check 3b | `avg_bound.py` |
| Tail form of the projection bound (§5.3); bootstrap intervals (§4.3, §5.1, §5.4, note to Table 3) | Check 3c | `bootstrap_tail_check.py` |
| First-order route share (§5.4; equation S2.1) | Exploratory | `route_share_check.py` |
| Exact check of the spread channel (Section S3) | Exploratory | `route_shock_check.py` |
| Table 4 (pickup tasks), §6 | Gate 4 | `lade_io.py`, `gate4.py` |
| §5.5 and Section S4 | Gate 5 | `gate5.py` (reuses the Gate 2c–4 results) |
| Figure 1 | — | `fig_exposure.py` |
| Theory: Lemmas 1–6, Propositions 1–7, Corollaries 1–4 (proofs in Section S1) | Checks T1–T18 | `theory_checks.py`, `test_sim.py` |
| Pre-specification record (Section S5) | — | `analysis/*/preregistration.md` and the specifications in `analysis/theory/` |

The figure scripts write a vector PDF and a 600 dpi PNG.

### Numerical checks in `theory_checks.py`

| Check | Result | Check | Result |
|---|---|---|---|
| T1 | Lemma 1 | T10 | Proposition 7 and Lemma 6 |
| T2 | Proposition 2 and Corollary 1 | T11 | Lemma 2 and the tail form of Proposition 5 |
| T3 | Proposition 3 (worst case at K = 1) | T12 | Proposition 2 attained |
| T4 | Proposition 5, route-level step | T13 | Propositions 1 and 3 on a bounded delay interval |
| T5, T6 | Proposition 4 and its tightness | T14 | Corollary 3 |
| T7 | Corollary 2 | T15 | Lemma 4 |
| T8 | Lemma 3 | T16 | Lemma 5 |
| T9 | Proposition 6 | T17 | Proposition 5 tight |
| | | T18 | Corollary 4 and its tightness |

All checks pass (`python -I code/resequencing/theory_checks.py`, about five minutes).

## Pre-specification

Each gate's `preregistration.md` was written and time-stamped before the gate's data were used. It logs any later deviation with its reason. These are project files, not entries in a public registry.

Two checks were specified after the gates:
- **Check 3b** (`analysis/theory/average_bound_check.md`), after Gate 3;
- **Check 3c and the bootstrap intervals** (`analysis/theory/projection_tail_and_bootstrap.md`), after Gate 5.

Each was specified before any of its quantities were computed.

## Notes on the derived data

In `analysis/gate4/gate4_routes.csv`, the LaDe courier identifiers are replaced by pseudonyms: their rank within the city, which preserves their order. The row order and all other columns are unchanged, so every result computed from the file is unchanged too.

The Amazon route identifiers and the Konovalenko et al. route identifiers are kept as published by the data providers.

## Licence

- **Code:** MIT (`LICENSE`).
- **Derived results and figures:** CC BY-NC 4.0 (`LICENSE-RESULTS.md`), following the licences of the source data.

Please cite the paper and this archive (`CITATION.cff`).
