# sbmlab-benchmarks

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23355028.svg)](https://doi.org/10.5281/zenodo.23355028)

Scripts, settings and result tables behind the measured posts on [sbmlab.com](https://sbmlab.com/blog). Each folder is one post: the code I ran, the outputs it printed, and a README with versions and the traps I hit.

Most proteomics folders use the same public three-species benchmark with known ratios, so every tool is scored against a right answer rather than against another tool.

Everything ran in October 2026 on one workstation: an i9-14900K with 128 GB RAM, Linux tools in WSL2 (64 GB available) and DIA-NN natively on Windows. Paths inside scripts and logs (`/home/sysoft/...`, `C:\diann-hye\...`) are from that machine; change them to yours.

## Posts and folders

### Proteomics: known-ratio benchmark (PRIDE PXD028735)

Human 65% in both conditions, yeast 30 → 15%, *E. coli* 5 → 20%, so the true log2(B/A) is 0 / −1 / +2. Three runs per condition.

| Folder | Post | What it shows |
|---|---|---|
| [maxquant-hye](maxquant-hye) | [MaxQuant 2.8 scored on known ratios](https://sbmlab.com/blog/maxquant-2-8-proteingroups-lfq-ibaq-mbr-scored-known-ratios-2026) | LFQ within ±0.5 of truth for 99 / 92 / 81% of human / yeast / *E. coli* proteins (raw intensity: 86 / 74 / 43%). Match-between-runs filled ~10% of peptide values; it made proteins quantified anyway more accurate, while proteins quantified only thanks to MBR were less reliable. The decoy column is now `Decoy`, not `Reverse` |
| [fragpipe-hye](fragpipe-hye) | [FragPipe vs MaxQuant](https://sbmlab.com/blog/fragpipe-vs-maxquant-2026-speed-benchmark) | Same 6 runs: 8.9 min by FragPipe's own timer (9.3 min wall clock in `time.log`, plus about 4 min raw → mzML in `convert.log`) vs 61.6 min. On the 3,877 proteins both quantified, within ±0.5: FragPipe 99.7 / 99.2 / 97.7% vs MaxQuant 98.9 / 92.2 / 80.9% |
| [diann-hye](diann-hye) | [DIA-NN on a known-answer benchmark](https://sbmlab.com/blog/diann-vs-maxquant-fragpipe-known-ratio-benchmark-2026) | DIA-NN 2.3.2 on the SWATH runs of the same samples (a different instrument): 8,987 protein groups; within ±0.5 98.6 / 95.1 / 85.4%. 6 h 17 min, ~3.5 h of it reading `.wiff` files |
| [limma-deqms-hye](limma-deqms-hye) | [limma vs DEqMS at small n](https://sbmlab.com/blog/limma-vs-deqms-proteomics-small-n-when-to-use-which-2026) | Unchanged human proteins called at adj. p < 0.05: t-test 68, limma 104, DEqMS 284 (an LFQ normalisation offset). With \|log2FC\| ≥ 0.5: 2 / 10 / 13 |
| [imputation-hye](imputation-hye) | [Imputation methods scored](https://sbmlab.com/blog/missing-values-proteomics-imputation-methods-compared) | Five imputation methods plus none. With a strict filter, imputation barely mattered. With a loose one it found more changes, doubled the false discoveries, and kNN erased a real 4-fold increase |
| [de-pipeline-check](de-pipeline-check) | [Differential expression in R](https://sbmlab.com/blog/differential-expression-proteomics-complete-r-pipeline) | Every R block of the post, run on the benchmark table (except the symbol-based enrichGO block) |
| [string-background](string-background) | [STRING enrichment background test](https://sbmlab.com/blog/string-enrichment-background-random-proteins-test-2026) | Random lists of detected proteins: median 103 "significant" GO/KEGG terms with STRING's default whole-genome background, 0 with the detected proteins as background |
| [r-vs-python-read-bench](r-vs-python-read-bench) | [R vs Python for bioinformatics](https://sbmlab.com/blog/r-vs-python-bioinformatics-guide-2026) | Reading and summarising a 146 MB MaxQuant `evidence.txt`: data.table 0.04 s, polars 0.16 s, pandas 0.51 s, base R 5.8 s |

### Preprint analyses

| Folder | What it adds |
|---|---|
| [preprint-hye](preprint-hye) | t-test / limma / DEqMS on the protein tables of all three engines (contaminants removed by accession), with and without an extra median normalisation; figure script for the preprint. Needs MaxQuant's `bin/conf/contaminants.fasta` copied next to the scripts as `contaminants_maxquant_2.8.1.fasta`; input paths at the top of each script point to the local MaxQuant, FragPipe and DIA-NN outputs |

### Other measured posts

| Folder | Post | Data | What it shows |
|---|---|---|---|
| [biomarker-leakage](biomarker-leakage) | [Data leakage with shuffled labels](https://sbmlab.com/blog/biomarker-discovery-data-leakage-shuffled-labels-auc-2026) | Bader et al. 2020 CSF proteome | Labels shuffled, proteins chosen before cross-validation: mean AUC 0.73. Done properly: 0.51 |
| [mrna-protein-wang2019](mrna-protein-wang2019) | [Proteomics vs genomics: mRNA vs protein](https://sbmlab.com/blog/proteomics-vs-genomics-differences) | Wang et al. 2019, 29 human tissues | mRNA explains ~29% of protein variation within a tissue; a gene's mRNA and protein track each other across tissues with median Spearman 0.35 |
| [bioinfo-db-check](bioinfo-db-check) | [13 bioinformatics databases, tested](https://sbmlab.com/blog/top-10-bioinformatics-databases-every-researcher-should-know) | Live public APIs | One protein (SOD1) through 13 APIs: what each returns, speed, failures, release, licence |
| [deseq2-edger-limma-airway](deseq2-edger-limma-airway) | [DESeq2 vs edgeR vs limma](https://sbmlab.com/blog/deseq2-vs-edger-vs-limma-rna-seq-comparison) | Bioconductor `airway` | The design (pairing by cell line) changed the results more than the choice of tool |
| [scrna-pbmc3k](scrna-pbmc3k) | [Single-cell RNA-seq on PBMC3k](https://sbmlab.com/blog/single-cell-rna-seq-analysis-complete-guide-2026) | 10x PBMC3k | Scanpy from QC to marker genes in ~23 s; two dependencies that are not installed by default |
| [biomart-pig-mouse](biomart-pig-mouse) | [BioMart pig ↔ mouse orthologues](https://sbmlab.com/blog/biomart-pig-mouse-ortholog-mapping-tutorial-2026) | Ensembl 111 and 116 | 1:1 orthologue tables, and why biomaRt failed where a plain martservice query worked |

A summary of all the tools, versions and timings is in [Best bioinformatics tools in 2026, installed and timed](https://sbmlab.com/blog/best-bioinformatics-tools-2026).

## What is not here

- **Raw data.** Download it from the source in each folder's README: PRIDE PXD028735, the journals' supplementary files, Bioconductor, 10x Genomics, Ensembl.
- **Third-party tables I only read**, such as the paper supplements and the raw API responses (KEGG does not allow redistribution). The scripts regenerate what they need.
- **Large logs and intermediate files.**

## Data sources and licences

- PXD028735: Van Puyvelde B. et al. (2022) *Scientific Data* 9:126, CC BY 4.0.
- Bader J.M. et al. (2020) *Molecular Systems Biology* 16:e9356, CC BY (open access, licence as listed by Europe PMC).
- Wang D. et al. (2019) *Molecular Systems Biology* 15:e8503, CC BY (open access, licence as listed by Europe PMC).
- `airway` (Himes et al. 2014) via Bioconductor; PBMC3k from 10x Genomics; Ensembl and STRING data under their own terms (see their websites).
- The STRING enrichment tables in `string-background` include KEGG pathway names returned by STRING.

Code in this repository is MIT licensed ([LICENSE](LICENSE)). Result tables derived from the datasets above keep those datasets' licences; please cite the original papers.

Found a mistake? Open an issue. These are single runs on one machine, so differences on your setup are expected and interesting.

## How to cite

Kang, B. (2026). *sbmlab-benchmarks: scripts and results behind the measured posts on sbmlab.com*. Zenodo. https://doi.org/10.5281/zenodo.23355028

This DOI always points to the latest version; v1.0.0 is https://doi.org/10.5281/zenodo.23355029. Please also cite the original datasets listed above.
