# Data (not redistributed)

Download the raw data into the folders below. Then build the processed files.

| Dataset | Source | Licence | Place the files in |
|---|---|---|---|
| 2021 Amazon Last Mile Routing Research Challenge (Merchán et al., 2024) | <https://registry.opendata.aws/amazon-last-mile-challenges> | CC BY-NC 4.0 | `data/raw/amazon_lmrrc_eval/` and `data/raw/amazon_lmrrc_train/` (the challenge's evaluation and training splits, with `route_data.json` in each) |
| Planned and driven routes (Konovalenko, Hvattum and Iversen, 2026) | <https://doi.org/10.17632/kkwgfvmtxn.1> | CC BY 4.0 | `data/raw/konovalenko2026/routes_performance.xlsx` |
| LaDe-P pickup data (Wu et al., 2024), five cities | <https://huggingface.co/datasets/Cainiao-AI/LaDe> | research use; CC BY-NC 4.0 per the dataset paper | `data/raw/lade/pickup/pickup_{cq,hz,jl,sh,yt}.csv` |

## Processed files

Run each command from the repository root. Each script states its exact arguments in its docstring.

| Output | Script |
|---|---|
| `data/processed/train_routes.pkl`, `data/processed/eval_routes.pkl` | `code/resequencing/amazon_io.py` |
| `data/processed/konovalenko_stops.pkl` | `python -I code/resequencing/konovalenko_io.py data/raw/konovalenko2026/routes_performance.xlsx data/processed/konovalenko_stops.pkl` |
| `data/processed/lade_pickup.pkl` | `code/resequencing/lade_io.py` |

`konovalenko_io.py` stores the spreadsheet as read. The eligibility rules and the split by week are applied downstream, as pre-specified.
