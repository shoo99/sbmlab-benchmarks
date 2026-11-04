#!/usr/bin/env bash
# Downloads the TCGA inputs of tcga_signatures.R into data/ next to this script (about 150 MB for the default cohorts).
#   <COHORT>_HiSeqV2.gz      Xena TCGA hub: gene-level RSEM, log2(norm + 1), Illumina HiSeq (UCSC Xena "TCGA.<COHORT>.sampleMap/HiSeqV2")
#   Survival_SupplementalTable_S1_20171025_xena_sp
#                            Xena PanCan Atlas hub: TCGA Pan-Cancer Clinical Data Resource (Liu et al. 2018, Cell 173:400),
#                            curated OS/DSS/DFI/PFI with AJCC stage, all 33 cancer types
# Usage: bash download.sh [COHORT ...]     (default: BRCA KIRC LUAD COAD)
# Files that already exist are kept; delete one to fetch it again. md5 sums go to data/md5.txt.
set -euo pipefail
cd "$(dirname "$0")" || exit 1
mkdir -p data
COHORTS=("$@"); [ ${#COHORTS[@]} -eq 0 ] && COHORTS=(BRCA KIRC LUAD COAD)
XENA_TCGA=https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download
XENA_PANCAN=https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download

fetch() {   # fetch <url> <file>
  if [ -s "data/$2" ]; then echo "have $2"; return; fi
  echo "get  $2"
  curl -fsSL --retry 3 -o "data/$2.part" "$1" && mv "data/$2.part" "data/$2"
}
fetch "$XENA_PANCAN/Survival_SupplementalTable_S1_20171025_xena_sp" Survival_SupplementalTable_S1_20171025_xena_sp
for c in "${COHORTS[@]}"; do
  fetch "$XENA_TCGA/TCGA.${c}.sampleMap%2FHiSeqV2.gz" "${c}_HiSeqV2.gz"
  gzip -t "data/${c}_HiSeqV2.gz"
done
(cd data && md5sum Survival_SupplementalTable_S1_20171025_xena_sp ./*_HiSeqV2.gz > md5.txt && cat md5.txt)
echo "next: Rscript tcga_signatures.R data ${COHORTS[*]} > tcga_signatures.out && Rscript score.R > score.out"
