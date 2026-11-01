# DIA-NN 2.7.0 vs 2.3.2 × PXD028735 HYE benchmark (2026-10-31 / 11-01)

Same six SCIEX TripleTOF 6600+ SWATH runs, same FASTA and same settings as [../diann-hye](../diann-hye) (DIA-NN 2.3.2), rerun with DIA-NN 2.7.0 Academia (Windows build, called from WSL, 24 threads, i9-14900K). Truth log2(B/A): human 0, yeast −1, E. coli +2.

Install: `DIA-NN-2.7.0-Academia.msi` from the GitHub release "DIA-NN 2.0", unpacked without installing (`msiexec /a … TARGETDIR=C:\DIA-NN\2.7.0`) so 2.3.2 stayed in place.

**.wiff trap:** the 2.7.0 package has no SCIEX reader DLLs. The first search (run_without_sciex_dlls.log) skipped all six .wiff files and still exited 0 after 2 min 35 s. Copying the 23 `*Clearcore*` / `*Sciex*` DLLs (originally from ProteoWizard, as the DIA-NN README says) from the 2.3.2 folder fixed it.

Run: `run270.sh` (library with the 2.3.2 parameters, 13 min 20 s; 6,370,129 precursors, 31,635 proteins, same as 2.3.2) and `search270.sh` (MBR search, `--temp` to a separate folder so the 2.3.2 .quant files were not overwritten; 9 h 10 min; the 12 .wiff loads took 326 min vs 199 min, +64%, other steps +25%; one run per version at different times, so version and machine cannot be separated). Logs: libgen.log, report.log.txt. Raw files are not included; download them from PRIDE.

Scoring (same rules as the 2.3.2 run: ≥2 of 3 per condition, mean log2, within ±0.5 of truth):
- `score270.R` → `score270.out` (report.pg_matrix.tsv and pg_mbr.tsv, the latter rebuilt from report.parquet with `../diann-hye/pg_from_parquet.py`)
- `common.R` → `common.out`: proteins quantified by both versions vs only by 2.7.0, with median intensity in condition A
- `filt.R` → `filt.out`: 2-of-3 vs 3-of-3 completeness rule, both versions
- `score270.R` also writes `diann_scores_270.csv`
- `export_fc.R` → `fc_both_versions.csv`, `chart270.py` → `diann-270-vs-232.png`

Main results:

| | 2.3.2 | 2.7.0 |
|---|---|---|
| Proteins per run (mean of 6, report.stats.tsv) | 8,198 | 8,708 |
| Quantified human / yeast / E. coli | 5,477 / 2,374 / 446 | 5,574 / 2,533 / 605 |
| Within ±0.5, all | 98.6 / 95.1 / 85.4% | 98.7 / 95.2 / 76.5% |
| Within ±0.5, the 8,097 shared | 99.0 / 95.9 / 86.1% | 99.2 / 97.2 / 87.2% |
| Within ±0.5, the 615 only in 2.7.0 | | 86.0 / 74.6 / 48.2% |

391 of the 615 groups quantified only by 2.7.0 were already in the 2.3.2 matrix with too few values (E. coli: 156 of 166, nearly all missing in condition A, where E. coli is 5%). They are low-intensity; their ratios scatter more and, for E. coli, are compressed toward 0 (median 1.48 instead of +2).

The default pg_matrix filter changed between the versions (2.5.1: no lib q-value filter with MBR; 2.6.0: 'Global' filters with MBR). With one filter for both (pg_mbr.tsv, Q ≤ 0.01 and PG.Q ≤ 0.01 from report.parquet), E. coli within ±0.5 goes from 87.6% of 421 (2.3.2) to 76.5% of 604 (2.7.0).

Scripts use paths relative to this folder and expect ../diann-hye next to it.

Post: https://sbmlab.com/blog/diann-2-7-0-upgrade-vs-2-3-2-known-ratio-benchmark-2026
