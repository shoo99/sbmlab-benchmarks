#!/usr/bin/env bash
# Downloads the inputs into data/ (not tracked; about 340 MB).
#   pbmc10k_multiome.h5  10x Genomics "PBMC from a healthy donor, granulocytes removed through cell sorting (10k)",
#                        Cell Ranger ARC 2.0.0 filtered feature-barcode matrix (RNA and ATAC peaks of the same nuclei),
#                        from 10x's public S3 bucket
#   hg38_genes.gtf       iGenomes UCSC hg38 gene annotation (for gene bodies and strands), from the iGenomes S3 bucket
# Usage: bash download.sh
set -euo pipefail
cd "$(dirname "$0")" || exit 1
mkdir -p data
TENX=https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k
IGEN=https://ngi-igenomes.s3.amazonaws.com/igenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes
[ -s data/pbmc10k_multiome.h5 ] || curl -fsSL --retry 3 -o data/pbmc10k_multiome.h5 "$TENX/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
[ -s data/hg38_genes.gtf ] || curl -fsSL --retry 3 -o data/hg38_genes.gtf "$IGEN/genes.gtf"
(cd data && md5sum pbmc10k_multiome.h5 hg38_genes.gtf | tee md5.txt)
