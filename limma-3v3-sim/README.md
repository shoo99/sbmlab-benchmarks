# limma-3v3-sim

Numbers behind [Three Samples per Group: What We Learned Switching from the t-Test to limma](https://sbmlab.com/blog/three-samples-per-group-t-test-vs-limma-omics).

`Rscript three_vs_three.R` (about 20 s) prints:

1. the smallest attainable two-sided Mann–Whitney p for 3, 4 and 5 samples per group (2 / C(2n, n));
2. the power of a two-sided two-sample t-test at n = 3 for d = 1, 2, 3, and the n per group needed for 80% power at d = 1 and 1.5 (`power.t.test`);
3. a synthetic experiment with a known answer: 400 features, 40 truly changed (log2 ±1.5, alternating up and down), 3 vs 3, per-feature within-group SD drawn from a log-normal with median 0.5 and sdlog 0.6 on the log2 scale, 500 repetitions, seed 1. Rules compared: p < 0.05; p < 0.05 and |log2FC| ≥ 1; BH-adjusted p < 0.05; each with Student's t-test and limma (`lmFit` + default `eBayes`).

`output.txt` is the printed output (R 4.3.3, limma 3.58.1). Synthetic data only; no real samples are involved.
