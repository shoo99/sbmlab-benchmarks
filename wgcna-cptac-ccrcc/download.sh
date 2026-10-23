#!/bin/bash
# Downloads the two public inputs (about 48 MB) into ./data
set -e
cd "$(dirname "$0")"; mkdir -p data/tjrcc
L=https://linkedomics.org/data_download/CPTAC-CCRCC
for f in HS_CPTAC_CCRCC_proteome_Tumor.cct HS_CPTAC_CCRCC_CLI.tsi; do curl -fsSL -o data/$f $L/$f; done
curl -fsSL -o data/tjrcc/matrix_protein.csv https://zenodo.org/api/records/8063124/files/matrix_protein.csv/content
md5sum -c md5sums.txt
