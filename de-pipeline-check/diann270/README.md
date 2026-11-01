# limma pipeline on DIA-NN 2.7.0 output, scored on known ratios (2026-11-01)

`diann_limma.R` runs the "Complete Minimal Pipeline" from https://sbmlab.com/blog/differential-expression-proteomics-complete-r-pipeline unchanged (except reading the input) on `../../diann-2.7-hye/report.pg_matrix.tsv` (PXD028735 SWATH, 3 vs 3, DIA-NN 2.7.0), then scores the calls: human true log2FC 0 (any call is false), yeast −1, E. coli +2. Output: `diann_limma.out`, `limma_results_filtered.csv`.

Main finding: the post's extra median normalization, applied on top of DIA-NN's own normalization, shifted every log2FC by about +0.18 in this design (more proteins go down than up), which pushed yeast to −0.88 and let only 339 of 2,533 yeast proteins pass |logFC| > 1. Without the extra step: 1,690 of 2,533, with fewer false human calls (3 vs 6).

With adj.P < 0.05 alone (no fold-change cutoff) the shift is worse: 3,320 of 6,218 calls (53.5%) are unchanged human proteins, against 173 of 3,066 (5.7%) without the extra normalization.
