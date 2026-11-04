#!/usr/bin/env bash
# Downloads Kang et al. 2018 (GEO GSE96583, batch 2: 8 donors, control and IFN-beta stimulated PBMC) into data/.
#   GSM2560248_2.1.mtx.gz, GSM2560248_barcodes.tsv.gz   control cells (kang_study.R uses only these)
#   GSM2560249_2.2.mtx.gz, GSM2560249_barcodes.tsv.gz   stimulated cells
#   GSE96583_batch2.genes.tsv.gz                         genes (Ensembl id, symbol), the rows of both matrices
#   GSE96583_batch2.total.tsne.df.tsv.gz                 per cell: donor (ind), condition (stim), cell type (cell),
#                                                        demuxlet call (multiplets: singlet / doublet / ambs)
# Usage: bash download.sh     (about 110 MB; GSE96583_RAW.tar is deleted after the four batch-2 files are extracted)
set -euo pipefail
cd "$(dirname "$0")" || exit 1
mkdir -p data
GEO=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl
fetch() {   # fetch <file>
  if [ -s "data/$1" ]; then echo "have $1"; return; fi
  echo "get  $1"
  curl -fsSL --retry 3 -o "data/$1.part" "$GEO/$1" && mv "data/$1.part" "data/$1"
}
fetch GSE96583_batch2.genes.tsv.gz
fetch GSE96583_batch2.total.tsne.df.tsv.gz
if [ ! -s data/GSM2560248_2.1.mtx.gz ] || [ ! -s data/GSM2560249_2.2.mtx.gz ]; then
  fetch GSE96583_RAW.tar
  tar -xf data/GSE96583_RAW.tar -C data --wildcards 'GSM256024[89]_*'
  rm -f data/GSE96583_RAW.tar
fi
for f in GSM2560248_2.1.mtx.gz GSM2560248_barcodes.tsv.gz GSM2560249_2.2.mtx.gz GSM2560249_barcodes.tsv.gz \
         GSE96583_batch2.genes.tsv.gz GSE96583_batch2.total.tsne.df.tsv.gz; do
  [ -s "data/$f" ] || { echo "missing data/$f"; exit 1; }
  gzip -t "data/$f"
done
(cd data && md5sum GSM256024[89]_* GSE96583_batch2.* > md5.txt && cat md5.txt)
echo "next: Rscript kang_study.R data > kang_study.out && Rscript score.R > score.out"
