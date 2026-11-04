# Unpaired RNA + ATAC integration scored with 10x Multiome pairing (2026-11-04)

Tracking issue: [#16](https://github.com/shoo99/sbmlab-benchmarks/issues/16). No post yet.

Methods that integrate single-cell RNA and ATAC from different cells are usually judged by how well the cell types mix. A 10x Multiome run measures both modalities in the same nuclei, so the right answer is known for every ATAC profile: it is the RNA profile with the same barcode. Here the pairing is hidden from each method and used only to score it.

TODO_GLUE_FINDINGS

## Data

- **Run**: 10x Genomics "PBMC from a healthy donor, granulocytes removed through cell sorting (10k)", Cell Ranger ARC 2.0.0 filtered feature-barcode matrix (`pbmc10k_multiome.h5`, md5 `df86844b99161b9487090d91e644745e`). It has 11,898 barcodes, 36,601 genes and 143,887 ATAC peaks.
- **Gene coordinates**: iGenomes UCSC hg38 `genes.gtf` (md5 `61535b7d94cdff8994d514e3b393cfe3`). Only transcripts that overlap the 10x gene's TSS interval are kept, so duplicated gene names do not create megabase spans.
- **Cells**: 11,310 pass the filters: RNA 1,000–25,000 counts, at least 500 genes, mitochondrial share < 20%, and ATAC 1,000–100,000 peak counts. 25 cells in an unassigned cluster are dropped, leaving **11,285**.
- **Downloads**: `download.sh` fetches both files. They are not tracked.

## Cell-type truth

The RNA side defines the labels. Each modality's own cells are never relabelled from the other.

- **Clustering**: log-normalised counts, 2,000 variable genes, 30 PCs, Louvain clusters at resolution 0.8.
- **Labels**: each cluster gets the type whose marker set has the highest mean scaled expression. A cluster whose best score is below 0.15 is left unassigned.

| Type | Markers | Nuclei |
|---|---|---|
| CD14 monocytes | CD14, LYZ, S100A8, VCAN | 3,180 |
| CD16 monocytes | FCGR3A, MS4A7, CDKN1C | 476 |
| cDC | FCER1A, CD1C, CLEC10A | 229 |
| pDC | LILRA4, IL3RA, CLEC4C | 111 |
| B | MS4A1, CD79A, BANK1, PAX5 | 962 |
| CD4 T | CD3E, IL7R, CD4, CCR7 | 3,188 |
| CD8 T | CD3E, CD8A, CD8B | 1,931 |
| NK | NKG7, GNLY, KLRD1, PRF1 | 1,208 |

The per-cluster scores are in `seurat_methods.out`.

## Methods

| Method | Version | What it does |
|---|---|---|
| Naive baseline | R 4.3.3, irlba 2.3.5.1 | ATAC gene activity (peak counts in the gene body plus 2 kb upstream, strand-aware) and RNA on the shared variable genes, each z-scored within its modality; one PCA (30 PCs) of both |
| Seurat CCA | Seurat 5.0.1 | `FindTransferAnchors(reduction = "cca")` from RNA to gene activity; `TransferData` for labels, weighted by the ATAC LSI (TF-IDF + SVD of peaks, dims 2–30); co-embedding as in Seurat's ATAC vignette (RNA imputed for the ATAC cells, merged, one PCA) |
| GLUE | scglue 0.5.2, torch TODO_TORCH (CPU) | the scglue tutorial: RNA counts with 2,000 Seurat-v3 variable genes and 100 PCs; ATAC LSI with 100 components; an RNA-anchored guidance graph from gene body + promoter to peaks; negative-binomial likelihoods; pretraining, then fine-tuning with the discriminator and balancing weights; random seed 0 |
| LIGER / UINMF | not run | `rliger` is not on conda, and CRAN is not reachable from the cloud environment |

- **Hidden pairing**: GLUE gets the ATAC cells under new names (`atac-0` …), so it cannot pair them by barcode. The R methods see the barcodes only as row names and never use them to match cells.
- **Labels**: every method transfers a type to each ATAC profile. Seurat uses `TransferData`. The baseline and GLUE use the majority of the 15 nearest RNA cells in their joint embedding.
- **Adding LIGER**: a run on a workstation can write `results/liger.tsv` in the same format (below), and `score.py` will pick it up.

## Scoring

Each `results/<method>.tsv` has one row per nucleus with these columns:
- `barcode`
- `true_type`: the RNA cell type
- `predicted_type`: the type transferred to its ATAC profile
- `foscttm_atac`: the share of RNA cells closer to the ATAC profile than its own RNA profile
- `foscttm_rna`: the reverse share

`score.py` (standard library) computes:

| Metric | Meaning |
|---|---|
| `foscttm` | fraction of samples closer than the true match, averaged over both directions (0 = every true match is the nearest; 0.5 = random) |
| `foscttm_atac_to_rna`, `foscttm_rna_to_atac` | each direction |
| `label_accuracy` | share of ATAC profiles given their own nucleus's RNA type |
| `errors_within_lineage` | share of the label errors that stay within a lineage (T/NK, B, monocytes, dendritic cells), e.g. CD4 T taken for CD8 T |
| per cell type | label accuracy and FOSCTTM |

## How to run

```bash
bash download.sh
```

```bash
Rscript seurat_methods.R > seurat_methods.out
```

```bash
python run_glue.py
```

```bash
python3 score.py > score.out && python3 figure.py
```

- **R methods**: `seurat_methods.R` took 7.4 minutes (4 BLAS threads) with a peak of 10.1 GB of memory. The 15 GB container ran out of memory twice before the script was changed to free the full matrix early and compute the LSI in place. It also writes `work/cells.tsv` and `work/gene_coords.tsv` for GLUE.
- **GLUE**: `run_glue.py` runs in a Python 3.11 environment with `pip install scglue scikit-misc igraph leidenalg` and pynndescent 0.5.13. The bioconda scglue pins numpy < 1.22, which cannot run scglue's own Python 3.10+ code.
  - With pynndescent 0.7.0, the balancing step fails: scanpy 1.11.5 calls `getnnz`, which the returned sparse array no longer has.
  - Each training stage is saved under `work/glue`, and a rerun resumes from the last saved stage.
  - The run here was resumed once, after pretraining, when the balancing step failed. Pretraining does not use pynndescent, igraph or leidenalg, the packages that changed in between. `run_glue.out` has both parts.
  - Run time: TODO_GLUE_TIME on 4 CPU threads.
- **CI**: re-runs `score.py` from the bundled `results/` and compares `scores.csv` and `score.out` byte for byte. It also runs `figure.py`.

## Results

![FOSCTTM and label-transfer accuracy per method and cell type](multiome_integration_pairing.png)

TODO_GLUE_RESULTS

## Traps

TODO_GLUE_TRAPS

## Files

| File | Content |
|---|---|
| `download.sh` | fetches the Multiome matrix and the gene annotation |
| `seurat_methods.R`, `seurat_methods.out` | cells, RNA cell types, gene activity, the naive baseline and Seurat CCA |
| `run_glue.py`, `run_glue.out` | GLUE |
| `results/` | one table per method: true type, transferred type and FOSCTTM per nucleus |
| `score.py`, `score.out`, `scores.csv` | scoring (standard library) |
| `figure.py`, `multiome_integration_pairing.png` | the figure |
