# Imputation methods on HYE known ratios (2026-10-10)
impute_bench.R <txt dir> <tag> → impute_<tag>.out / .csv.gz. MBR on: ~/mq-hye/results_mbr_on/txt; MBR off: ~/mq-hye/raw/combined/txt.
chart.py → sbmlab public/images/blog/imputation-hye-ecoli.png. msImpute install failed on R 4.3 (not tested).

## How it scores

- Input: MaxQuant 2.8.1 `LFQ intensity`, decoys, contaminants and site-only groups removed, species from the `OX=` taxon IDs of all proteins in the group (one species only).
- `E.coli median lfc`, `yeast median lfc` and `within ±0.5` use log2FC **centred on the human median** (`lfc_c`). This differs from `maxquant-hye/score_quant.R`, which scores uncentred log2FC. On the same 2,896 / 907 / 132 proteins (`none`, MBR on), yeast within ±0.5 is 97% centred here vs 92.3% uncentred in `maxquant-hye/quant_scores.csv`. Don't put the two side by side without saying which rule each uses.
- A protein is called at limma adj.P < 0.05 and |logFC| ≥ 0.5, where logFC is limma's uncentred log2(B/A). A call counts as found only when it goes the right way (yeast down, E. coli up). Human calls and wrong-direction calls are both false, and the FDP is (human + wrong direction) / calls.

## Correction (2026-11-02): wrong-direction calls

The version that produced the bundled `impute_<tag>.out` and `.csv.gz` counted every non-human call as found, including yeast called up and E. coli called down, and left those calls out of the FDP. `impute_bench.R` now scores by direction and also writes `impute_<tag>_summary.csv`. The raw MaxQuant tables are not in the repository, so the bundled outputs were not re-run. Instead, `direction_rescore.py` re-scores them: it recovers each method's human-median offset from the printed counts, and checks that the result is unique and that the `none` offset matches `maxquant-hye/quant_scores.csv` (+0.144). It writes `impute_<tag>_summary.csv` and `direction_rescore.out`.

With the strict filter (≥2 of 3 in both groups), no method had a wrong-direction call, so those rows are unchanged. With the loose filter (≥2 of 3 in either group):

| Data | Method | Found (before → after) | FDP (before → after) |
|---|---|---|---|
| MBR on | MinDet | 1102 → 1100 of 1415 | 2.3% → 2.5% |
| MBR on | MinProb | 1098 → 1096 | 2.7% → 2.8% |
| MBR on | kNN | 977 → **961** | 2.9% → **4.5%** (16 wrong-direction calls) |
| MBR off | MinDet / MinProb / Perseus | −1 each | +0.1 points each |
| MBR off | kNN | 778 → **755** of 1163 | 3.7% → **6.6%** (23 wrong-direction calls) |

Perseus and QRILC (MBR on) and QRILC (MBR off) are unchanged. The conclusion holds and gets stronger: with the loose filter, imputation finds more changes but at 2–5 times the strict-filter FDP of the same method, and kNN is the worst. `preprint-hye/figures.py` (Fig 3a–b) now reads the summary files.
