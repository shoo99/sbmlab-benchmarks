# mRNA vs protein in 29 human tissues (2026-10-14)
Data: Wang D. et al. (2019) Mol Syst Biol 15:e8503, Table EV6 sheet "A. Protein vs mRNA correlation"
(log10 FPKM and log10 iBAQ, 9,870 genes detected as both in >= 10 tissues). The xlsx is not included; download it from the journal's supplementary files.

- `analyse.py Table_EV6.xlsx` → `analyse.out`, `within_tissue.csv`, `per_gene.csv` (all numbers in the post).
- `mrna_protein.py` = the code block in the sbmlab post, run as written (output matches the post).
- `chart.py` → the post's figure.

Headline: within a tissue, mRNA vs protein Spearman median 0.52 (0.41–0.58), R² median 0.29;
per gene across tissues, median Spearman 0.35. Per-gene values agree with the paper's own column (Pearson 1.000).
