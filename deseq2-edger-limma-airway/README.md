# DESeq2 vs edgeR (QL) vs limma-voom on airway (2026-10-10)
`airway_three.R` → `airway_three.out`. R 4.3.3; DESeq2 1.42.1, edgeR 4.0.16, limma 3.58.1, airway 1.22.0.
Same code as the sbmlab post (minus lfcShrink), run with ~condition (as written) and ~cell + condition (paired).
FDR<0.05: ~condition 2,782/1,952/1,878 (all three 1,811, 65% of union); ~cell+condition 4,099/5,130/5,170 (4,055, 77%).
|log2FC|≥1: 828/661/635 (594, 69%) and 1,000/1,169/1,155 (960, 79%). Pairing: DESeq2 +1,414 / −97. Each tool 0.4–3.1 s.
