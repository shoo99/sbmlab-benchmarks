# limma pipeline on DIA-NN 2.7.0 output, scored on known ratios (2026-11-01)

**2026-11-23: re-run on verified raw files.** The first PRIDE downloads of 5 of the 6 SWATH `.wiff.scan` files contained short server error messages (14 damaged spectra). DIA-NN 2.7.0 was re-run on verified copies with the same settings, `diann_limma.R` was re-run on the new `report.pg_matrix.tsv`, and every number below is from the re-run. Old → new: proteins after the per-group filter 8,712 → 8,506; yeast passing |logFC| > 1 with the extra normalization 339 of 2,533 → 604 of 2,466, without it 1,690 → 1,727; false human share with adj.P alone 53.4% → 51.3%.

`diann_limma.R` runs the "Complete Minimal Pipeline" from https://sbmlab.com/blog/differential-expression-proteomics-complete-r-pipeline unchanged (except reading the input) on `../../diann-2.7-hye/report.pg_matrix.tsv` (PXD028735 SWATH, 3 vs 3, DIA-NN 2.7.0), then scores the calls: human true log2FC 0 (any call is false), yeast −1, E. coli +2. It then re-runs the same steps with variations (below). Output: `diann_limma.out`, `limma_results_filtered.csv`.

The pipeline exactly as written, with no missing-value filter, tests 8,727 proteins and calls 1,149: 562 of 2,552 yeast down, 573 of 614 E. coli up, 12 human (`diann_limma.out` line 4). The numbers below add the post's Step 1 per-group filter (≥ 2 of 3 values in each group) before normalising, which leaves 8,506 proteins (lines 6–14).

Main finding: the post's extra median normalization, applied on top of DIA-NN's own normalization, shifted every log2FC by about +0.18 in this design (more proteins go down than up). That pushed yeast to −0.91 and let only 604 of 2,466 yeast proteins pass |logFC| > 1. Without the extra step: 1,727 of 2,466, with fewer false human calls (5 vs 7).

With adj.P < 0.05 alone (no fold-change cutoff) the shift is worse: 2,898 of 5,651 calls (51.3%) are unchanged human proteins, against 159 of 2,899 (5.5%) without the extra normalization. The `false share` printed in `diann_limma.out` (51.4% and 5.6%) also counts the 5 and 2 calls in the wrong direction.

2026-11-02: species assignment now follows the rule used in the other folders, where a group naming two species is left unassigned instead of being called Human. This matrix has no such group. A re-run with R 4.3.3 and limma 3.58.1 reproduced `diann_limma.out` and `limma_results_filtered.csv` byte for byte, before and after the change.
