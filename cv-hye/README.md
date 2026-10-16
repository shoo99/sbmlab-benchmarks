# cv-hye — replicate CV on the HYE benchmark, three engines

Post: [Proteomics CV: what replicate CVs look like in MaxQuant, FragPipe and DIA-NN](https://sbmlab.com/blog/proteomics-cv-replicates-maxquant-fragpipe-diann-measured-2026)

Data: PRIDE PXD028735 (Van Puyvelde et al. 2022), conditions A and B, preparations Alpha/Beta/Gamma, injection 01 each.
Engine outputs are the ones from the `maxquant-hye`, `fragpipe-hye` and `diann-hye` folders (MBR on). Input paths at the top of `cv_hye.R` point to the local result files; change them to yours.

| File | What it does |
|---|---|
| `cv_hye.R` → `cv_hye.out`, `cv_summary.csv`, `cv_by_quartile.csv`, `cv_vs_accuracy.csv` | CV = sd/mean of the three linear intensities per protein and condition (only when all three are present); MaxQuant LFQ vs raw Intensity; CV by abundance quartile (within engine, condition A); ratio accuracy (±0.5 log2 of truth after centring on the human median) split by the larger of the two CVs at 20% |
| `completeness_check.R` → `completeness_check.out` | Share of MaxQuant proteins with all three LFQ values, by raw-intensity quartile (condition A) |
| `figure.py` → `cv_by_quartile.png` | Figure in the post |

Contaminants: groups containing any accession from MaxQuant 2.8.1 `bin/conf/contaminants.fasta` are removed from FragPipe and DIA-NN (copy it next to the script as `contaminants_maxquant_2.8.1.fasta`).
The replicates are preparation + injection replicates, not biological replicates.
