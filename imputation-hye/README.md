# Imputation methods on HYE known ratios (2026-10-10)
impute_bench.R <txt dir> <tag> → impute_<tag>.out / .csv.gz. MBR on: ~/mq-hye/results_mbr_on/txt; MBR off: ~/mq-hye/raw/combined/txt.
chart.py → sbmlab public/images/blog/imputation-hye-ecoli.png. msImpute install failed on R 4.3 (not tested).

**2026-11-23: re-run on verified raw files.** The downloader had written a ~200-byte server error message into one spot each of A_Alpha_01 and A_Gamma_01. Clean copies were downloaded and MaxQuant was re-run with the same settings (MBR on and off); all `impute_<tag>` outputs are from the new runs. MBR on, strict filter: 3,935 → 3,941 proteins, `none` 958 of 1,039 → 956 of 1,045 found. Loose-filter kNN FDP: 4.5% → 4.3% (MBR on), 6.6% → 5.0% (MBR off). MBR off moved more because the old MBR-off run did not match the MBR-on run in per-run peptide counts; the new runs match.

## How it scores

- Input: MaxQuant 2.8.1 `LFQ intensity`, decoys, contaminants and site-only groups removed, species from the `OX=` taxon IDs of all proteins in the group (one species only).
- `E.coli median lfc`, `yeast median lfc` and `within ±0.5` use log2FC **centred on the human median** (`lfc_c`). This differs from `maxquant-hye/score_quant.R`, which scores uncentred log2FC. On the same 2,896 / 912 / 133 proteins (`none`, MBR on), yeast within ±0.5 is 97% centred here vs 91.6% uncentred in `maxquant-hye/quant_scores.csv`. Don't put the two side by side without saying which rule each uses.
- A protein is called at limma adj.P < 0.05 and |logFC| ≥ 0.5, where logFC is limma's uncentred log2(B/A). A call counts as found only when it goes the right way (yeast down, E. coli up). Human calls and wrong-direction calls are both false, and the FDP is (human + wrong direction) / calls.

## Correction (2026-11-02): wrong-direction calls

The version that produced the bundled `impute_<tag>.out` and `.csv.gz` counted every non-human call as found, including yeast called up and E. coli called down, and left those calls out of the FDP. `impute_bench.R` now scores by direction and also writes `impute_<tag>_summary.csv`. The raw MaxQuant tables are not in the repository, so the bundled outputs were not re-run. Instead, `direction_rescore.py` re-scores them: it recovers each method's human-median offset from the printed counts, and checks that the result is unique and that the `none` offset matches `maxquant-hye/quant_scores.csv` (+0.144). It writes `impute_<tag>_summary.csv` and `direction_rescore.out`.

**Re-run on the raw data (2026-11-03).** The fixed `impute_bench.R` was then re-run on the MaxQuant 2.8.1 tables (`~/mq-hye/...`). All 24 summary rows came out identical to the re-scored ones, and the centred log2FC and adj.P of all 50,215 / 40,768 rows match the previous `.csv.gz`. The bundled `impute_<tag>.out`, `.csv.gz` (now with a `logFC` column) and `_summary.csv` were from that re-run until they were replaced by the 2026-11-23 re-run (note at the top). `direction_rescore.py` and `direction_rescore.out` are kept as the record of the re-score; the script reads the outputs of the old version, which are in the git history before this change.

The table below and `direction_rescore.out` predate the 2026-11-23 re-run and describe the old search. With the strict filter (≥2 of 3 in both groups), no method had a wrong-direction call, so those rows are unchanged. With the loose filter (≥2 of 3 in either group):

| Data | Method | Found (before → after) | FDP (before → after) |
|---|---|---|---|
| MBR on | MinDet | 1102 → 1100 of 1415 | 2.3% → 2.5% |
| MBR on | MinProb | 1098 → 1096 | 2.7% → 2.8% |
| MBR on | kNN | 977 → **961** | 2.9% → **4.5%** (16 wrong-direction calls) |
| MBR off | MinDet / MinProb / Perseus | −1 each | +0.1 points each |
| MBR off | kNN | 778 → **755** of 1163 | 3.7% → **6.6%** (23 wrong-direction calls) |

Perseus and QRILC (MBR on) and QRILC (MBR off) are unchanged.

Current values (2026-11-23 re-run, `impute_<tag>_summary.csv`). The strict filter still has no wrong-direction calls. Loose filter:

| Data | Method | Found | FDP |
|---|---|---|---|
| MBR on | MinDet | 1092 of 1417 | 2.5% |
| MBR on | MinProb | 1089 | 2.5% |
| MBR on | Perseus | 1038 | 2.4% |
| MBR on | QRILC | 1017 | 2.3% |
| MBR on | kNN | 996 | **4.3%** (19 wrong-direction calls) |
| MBR off | MinDet | 867 of 1183 | 3.1% |
| MBR off | MinProb | 868 | 3.2% |
| MBR off | Perseus | 818 | 3.5% |
| MBR off | QRILC | 783 | 3.6% |
| MBR off | kNN | 743 | **5.0%** (15 wrong-direction calls) |

The conclusion holds: with the loose filter, imputation finds more changes but at 2 to 7 times the strict-filter FDP of the same method (2 to 5 times before the re-run; the top of the range is QRILC without MBR, 0.5% → 3.6%), and kNN is the worst. `preprint-hye/figures.py` (Fig 3a–b) now reads the summary files.
