#!/usr/bin/env bash
# Downloads the two doublet truth sets into data/ (about 220 MB).
#   Kang et al. 2018, GEO GSE96583 batch 2 (demuxlet doublet calls between 8 pooled donors):
#     GSM2560248_2.1.mtx.gz, GSM2560248_barcodes.tsv.gz   control lane
#     GSM2560249_2.2.mtx.gz, GSM2560249_barcodes.tsv.gz   IFN-beta lane
#     GSE96583_batch2.genes.tsv.gz, GSE96583_batch2.total.tsne.df.tsv.gz   genes; per-cell donor, lane, demuxlet call
#   10x Genomics hgmm_12k (Cell Ranger 2.1.0): about 12,000 droplets of mixed HEK293T (human) and NIH3T3 (mouse) cells
#     hgmm_12k/filtered_gene_bc_matrices/{hg19,mm10}/
# Usage: bash download.sh      (files that already exist are kept)
set -euo pipefail
cd "$(dirname "$0")" || exit 1
mkdir -p data
GEO=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl
TENX=https://cf.10xgenomics.com/samples/cell-exp/2.1.0/hgmm_12k/hgmm_12k_filtered_gene_bc_matrices.tar.gz
fetch() {   # fetch <url> <file>
  if [ -s "data/$2" ]; then echo "have $2"; return; fi
  echo "get  $2"
  curl -fsSL --retry 3 -o "data/$2.part" "$1" && mv "data/$2.part" "data/$2"
}
fetch "$GEO/GSE96583_batch2.genes.tsv.gz" GSE96583_batch2.genes.tsv.gz
fetch "$GEO/GSE96583_batch2.total.tsne.df.tsv.gz" GSE96583_batch2.total.tsne.df.tsv.gz
if [ ! -s data/GSM2560248_2.1.mtx.gz ] || [ ! -s data/GSM2560249_2.2.mtx.gz ]; then
  fetch "$GEO/GSE96583_RAW.tar" GSE96583_RAW.tar
  tar -xf data/GSE96583_RAW.tar -C data --wildcards 'GSM256024[89]_*'
  rm -f data/GSE96583_RAW.tar
fi
if [ ! -d data/hgmm_12k/filtered_gene_bc_matrices ]; then
  fetch "$TENX" hgmm_12k_filtered_gene_bc_matrices.tar.gz
  mkdir -p data/hgmm_12k
  tar -xzf data/hgmm_12k_filtered_gene_bc_matrices.tar.gz -C data/hgmm_12k
  rm -f data/hgmm_12k_filtered_gene_bc_matrices.tar.gz
fi
ls data/hgmm_12k/filtered_gene_bc_matrices/hg19/ data/hgmm_12k/filtered_gene_bc_matrices/mm10/ > /dev/null
(cd data && md5sum GSM256024[89]_* GSE96583_batch2.* hgmm_12k/filtered_gene_bc_matrices/*/* > md5.txt && cat md5.txt)
echo "next: bash run.sh    (prepare, run every tool, collect and score)"
