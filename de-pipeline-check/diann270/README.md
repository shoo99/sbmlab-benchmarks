# limma pipeline on DIA-NN 2.7.0 output, scored on known ratios (2026-11-01)

`diann_limma.R` runs the "Complete Minimal Pipeline" from https://sbmlab.com/blog/differential-expression-proteomics-complete-r-pipeline unchanged (except reading the input) on `../../diann-2.7-hye/report.pg_matrix.tsv` (PXD028735 SWATH, 3 vs 3, DIA-NN 2.7.0), then scores the calls: human true log2FC 0 (any call is false), yeast −1, E. coli +2. It then re-runs the same steps with variations (below). Output: `diann_limma.out`, `limma_results_filtered.csv`.

The pipeline exactly as written, with no missing-value filter, tests 8,905 proteins and calls 1,002: 392 of 2,622 yeast down, 596 of 643 E. coli up, 10 human (`diann_limma.out` line 4). The numbers below add the post's Step 1 per-group filter (≥ 2 of 3 values in each group) before normalising, which leaves 8,712 proteins (lines 6–14).

Main finding: the post's extra median normalization, applied on top of DIA-NN's own normalization, shifted every log2FC by about +0.18 in this design (more proteins go down than up). That pushed yeast to −0.88 and let only 339 of 2,533 yeast proteins pass |logFC| > 1. Without the extra step: 1,690 of 2,533, with fewer false human calls (3 vs 6).

With adj.P < 0.05 alone (no fold-change cutoff) the shift is worse: 3,320 of 6,218 calls (53.4%) are unchanged human proteins, against 173 of 3,066 (5.6%) without the extra normalization. The `false share` printed in `diann_limma.out` (53.5% and 5.7%) also counts the 8 and 2 calls in the wrong direction.

2026-11-02: species assignment now follows the rule used in the other folders, where a group naming two species is left unassigned instead of being called Human. This matrix has no such group. A re-run with R 4.3.3 and limma 3.58.1 reproduced `diann_limma.out` and `limma_results_filtered.csv` byte for byte, before and after the change.
