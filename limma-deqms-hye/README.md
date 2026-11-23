# t-test vs limma vs DEqMS on the HYE known-ratio benchmark (2026-10-10)
Input: MaxQuant 2.8.1 proteinGroups.txt with MBR, PXD028735 Orbitrap DDA 3 vs 3. The current outputs come from the re-run table (G:\dda-clean\mq\mq_target_only_mbr_on); the path in the script, ~/mq-hye/results_mbr_on/txt, still holds the first run.
`limma_deqms_hye.R` → `limma_deqms_hye.out`, `limma_deqms_hye_results.csv`. R 4.3.3, limma, DEqMS 1.20.0.

**2026-11-23: re-run on verified raw files.** The first PRIDE downloads of A_Alpha_01 and A_Gamma_01 each contained a ~200-byte server error message (one damaged spectrum per file). MaxQuant was re-run on verified copies with the same settings, this script was re-run on the new proteinGroups.txt, and every number below is from the re-run. Old → new at adj.p<0.05: proteins tested 3,532 → 3,531, limma true/false calls 891/104 → 889/101, DEqMS 893/284 → 891/275.

3,531 proteins complete in both groups (2,626 human / 797 yeast / 108 E. coli); only 8 single-peptide.
adj.p<0.05: t 508/905 found, 67 human FP · limma 889, 101 · DEqMS 891, 275 (99% of DEqMS human calls positive; human median log2FC +0.144 = LFQ normalisation offset).
+|log2FC|≥0.5: t 504/2 · limma 873/10 · DEqMS 877/13.
`post_code_check.R`: the sbmlab post's code blocks extracted verbatim (only file path + 6-sample group line changed) → runs; 5,628 rows after filters.

2026-11-02: lines 17–20 of `limma_deqms_hye.out` (human median log2FC, and per method the human calls' median |lfc| and share with lfc > 0) were printed by code missing from the committed script. That code is now back at the end of `limma_deqms_hye.R`. Run on `limma_deqms_hye_results.csv`, it reproduces those lines exactly. `score()` now also requires the right direction for a true call (yeast down, E. coli up). No call here went the wrong way, so all printed numbers are unchanged; the scoring part, run on the results CSV, reproduces lines 5–20.
