# DIA-NN quantification modes on known ratios

Blog post: [DIA-NN quantification modes scored on known ratios](https://sbmlab.com/blog/diann-quantums-vs-legacy-maxlfq-fragment-quantification-known-ratios-2026)

**2026-11-23: re-run on verified raw files.** The first download had written ~200-byte server error messages into 5 of the 6 .wiff.scan files (one spectrum per spot). The six DIA-NN runs, `logs/c_*.report.log.txt` and everything in `clean/` come from a re-run on clean copies (directLFQ 0.3.3, iq 2.0.1). The top-level files were not re-run.
*E. coli* within ±0.5, MBR off: QuantUMS default 89.6% → 91.9%, legacy MaxLFQ 91.6% → 93.3%, fragment directLFQ 90.8% → 89.1%. Without MBR, directLFQ now scores below both QuantUMS modes instead of above them.

PXD028735 (Van Puyvelde et al. 2022), TripleTOF 6600+ SWATH, conditions A and B, first injection of Alpha/Beta/Gamma (3 vs 3).
Truth log2(B/A): human 0, yeast −1, *E. coli* +2. DIA-NN 2.3.2, library predicted from the FASTA, as in [diann-hye](../diann-hye).

## Runs

`logs/clean_rerun.sh` (log: `logs/clean_rerun.log`; DIA-NN logs: `logs/c_*.report.log.txt`). Mass accuracies and scan window fixed
(`--mass-acc 16 --mass-acc-ms1 13 --window 9`), `--export-quant` for fragment intensities.
The "exit 0" at the end of each run in `logs/clean_rerun.log` is meaningless: the script read `$?` after `$(date)` had reset it (fixed 2026-11-02). Each run's own log, `logs/c_*.report.log.txt`, ends with "Stats report saved".
`logs/clean_rerun.log` is from the original October run on the damaged raw files. The 2026-11-23 re-run used the same DIA-NN commands with only the raw folder, output folders and `--temp` changed; its DIA-NN logs replaced `logs/c_*.report.log.txt`, and the wall times below are the last timestamps in those logs.

| Run | Flags | Wall time |
|---|---|---|
| c_mbr_quantums | `--reanalyse` (writes the first-pass `.quant` files) | 9 h 0 min |
| c_nombr_quantums / highacc / direct | `--use-quant` (+ `--high-acc` / `--direct-quant`) | 4 / 4 / 2 min |
| c_mbr_highacc / direct | `--use-quant --reanalyse` (+ `--high-acc` / `--direct-quant`) | 2 h 44 min / 2 h 41 min |

## Scripts (run in this order)

| Script | What it does |
|---|---|
| `build_inputs.py <set>` | `set` = `nombr`, `mbr`, `nombr_runq` (or `first`). Writes `clean/<set>/`: DIA-NN's `pg_matrix` per mode, precursor sums, and `fragments_long.tsv` (fragment apex intensities × DIA-NN's per-precursor normalisation factor). Precursor filter reproduces DIA-NN's own `pr_matrix` exactly (nombr: Global q; mbr: Lib q; nombr_runq: run-specific Q.Value/PG.Q.Value, kept to show the filter effect) |
| `fragment_maxlfq_iq.R clean/<set>` | MaxLFQ on fragments with iq (median normalisation off) |
| `fragment_directlfq.py clean/<set>` | directLFQ on the same fragments |
| `score_modes.R clean/<set>` | Scores every route: ≥ 2 of 3 values per condition, within ±0.5 of truth, IQR, CV. Output `score_modes.out`, `scores_all.csv`, `scores_common.csv` |
| `mbr_split.R` | MBR off vs on, split by proteins also quantified without MBR / only with MBR; MBR-added *E. coli* split by precursor completeness. Output `clean/mbr_split.out` and CSVs |
| `figure.py` | `clean/ecoli_routes.png` |

`score_modes.R`, `mbr_split.R` and `figure.py` also need MaxQuant 2.8.1's `bin/conf/contaminants.fasta` (path at the top of each script; not in this repository). They drop every protein group that contains one of its accessions, matched exactly as written in the FASTA header, so an entry with an isoform suffix (`-1`) does not match the canonical accession DIA-NN reports. Without contaminant removal, clean/nombr scores 4 more human proteins (QuantUMS default: 5,099 vs 5,095); *E. coli* is unchanged.

The fragment tables (`fragments_long.tsv`, `fragments_wide.aq_reformat.tsv`, `*.ion_intensities.tsv`; hundreds of MB) are not committed; `build_inputs.py` recreates them from DIA-NN's `report.parquet`.

Files in the top level of this folder (`pg_*.tsv`, `score_modes.out`, `scores_*.csv`, `directlfq.log`, `fragments_wide*`) are from a first exploratory rerun with automatic mass accuracy and the run-specific filter (`build_inputs.py first`). That run used the damaged raw files and was not repeated on 2026-11-23, so these files are kept as they were. The blog post uses the `clean/` results. They were built before `build_inputs.py` started writing `pg_precsum6_*.tsv`, so the top-level `score_modes.out` has no "precursors in all 6 runs" row. Re-running `build_inputs.py first` and `score_modes.R` would add it.

## Main results (*E. coli*, within ±0.5 of truth)

| Route | MBR off | MBR on, also without MBR | MBR on, only with MBR |
|---|---|---|---|
| QuantUMS default | 91.9% | 92.9% | 72.2% |
| QuantUMS high-accuracy | 91.9% | 93.6% | 72.2% |
| Legacy MaxLFQ | 93.3% | 96.3% | 62.8% |
| Fragment directLFQ | 89.1% | 93.0% | 48.6% |
| Fragment MaxLFQ (iq) | 89.1% | 86.4% | 40.7% |
| Precursor sum | 33.1% | 51.5% | 25.0% |
