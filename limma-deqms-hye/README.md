# t-test vs limma vs DEqMS on the HYE known-ratio benchmark (2026-10-10)
Input: MaxQuant 2.8.1 proteinGroups.txt with MBR (~/mq-hye/results_mbr_on/txt), PXD028735 Orbitrap DDA 3 vs 3.
`limma_deqms_hye.R` → `limma_deqms_hye.out`, `limma_deqms_hye_results.csv`. R 4.3.3, limma, DEqMS 1.20.0.
3,532 proteins complete in both groups (2,626 human / 798 yeast / 108 E. coli); only 7 single-peptide.
adj.p<0.05: t 507/906 found, 68 human FP · limma 891, 104 · DEqMS 893, 284 (99% of DEqMS human calls positive; human median log2FC +0.145 = LFQ normalisation offset).
+|log2FC|≥0.5: t 503/2 · limma 877/10 · DEqMS 880/13.
`post_code_check.R`: the sbmlab post's code blocks extracted verbatim (only file path + 6-sample group line changed) → runs; 5,628 rows after filters.
