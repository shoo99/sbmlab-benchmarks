# DIA-NN quantification modes on known ratios

Blog post: [DIA-NN quantification modes scored on known ratios](https://sbmlab.com/blog/diann-quantums-vs-legacy-maxlfq-fragment-quantification-known-ratios-2026)

PXD028735 (Van Puyvelde et al. 2022), TripleTOF 6600+ SWATH, conditions A and B, first injection of Alpha/Beta/Gamma (3 vs 3).
Truth log2(B/A): human 0, yeast −1, *E. coli* +2. DIA-NN 2.3.2, library predicted from the FASTA, as in [diann-hye](../diann-hye).

## Runs

`logs/clean_rerun.sh` (log: `logs/clean_rerun.log`; DIA-NN logs: `logs/c_*.report.log.txt`). Mass accuracies and scan window fixed
(`--mass-acc 16 --mass-acc-ms1 13 --window 9`), `--export-quant` for fragment intensities.

| Run | Flags | Wall time |
|---|---|---|
| c_mbr_quantums | `--reanalyse` (writes the first-pass `.quant` files) | 9 h 20 min |
| c_nombr_quantums / highacc / direct | `--use-quant` (+ `--high-acc` / `--direct-quant`) | 4 / 4 / 2 min |
| c_mbr_highacc / direct | `--use-quant --reanalyse` (+ `--high-acc` / `--direct-quant`) | 2 h 57 min / 3 h 6 min |

## Scripts (run in this order)

| Script | What it does |
|---|---|
| `build_inputs.py <set>` | `set` = `nombr`, `mbr`, `nombr_runq` (or `first`). Writes `clean/<set>/`: DIA-NN's `pg_matrix` per mode, precursor sums, and `fragments_long.tsv` (fragment apex intensities × DIA-NN's per-precursor normalisation factor). Precursor filter reproduces DIA-NN's own `pr_matrix` exactly (nombr: Global q; mbr: Lib q; nombr_runq: run-specific Q.Value/PG.Q.Value, kept to show the filter effect) |
| `fragment_maxlfq_iq.R clean/<set>` | MaxLFQ on fragments with iq (median normalisation off) |
| `fragment_directlfq.py clean/<set>` | directLFQ on the same fragments |
| `score_modes.R clean/<set>` | Scores every route: ≥ 2 of 3 values per condition, within ±0.5 of truth, IQR, CV. Output `score_modes.out`, `scores_all.csv`, `scores_common.csv` |
| `mbr_split.R` | MBR off vs on, split by proteins also quantified without MBR / only with MBR; MBR-added *E. coli* split by precursor completeness. Output `clean/mbr_split.out` and CSVs |
| `figure.py` | `clean/ecoli_routes.png` |

The fragment tables (`fragments_long.tsv`, `fragments_wide.aq_reformat.tsv`, `*.ion_intensities.tsv`; hundreds of MB) are not committed; `build_inputs.py` recreates them from DIA-NN's `report.parquet`.

Files in the top level of this folder (`pg_*.tsv`, `score_modes.out`, `scores_*.csv`, `directlfq.log`) are from a first exploratory rerun with automatic mass accuracy and the run-specific filter (`build_inputs.py first`). The blog post uses the `clean/` results.

## Main results (*E. coli*, within ±0.5 of truth)

| Route | MBR off | MBR on, also without MBR | MBR on, only with MBR |
|---|---|---|---|
| QuantUMS default | 89.6% | 94.6% | 74.0% |
| QuantUMS high-accuracy | 90.3% | 95.3% | 70.5% |
| Legacy MaxLFQ | 91.6% | 96.6% | 62.4% |
| Fragment directLFQ | 90.8% | 93.0% | 45.1% |
| Fragment MaxLFQ (iq) | 87.5% | 87.4% | 33.1% |
| Precursor sum | 34.3% | 53.6% | 25.7% |
