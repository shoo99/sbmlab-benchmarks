#!/usr/bin/env bash
# Downloads for go-release-drift (issue #19); nothing here is tracked. About 200 MB.
#   GO releases (ontology go-basic.obo + human GAF) from the Gene Ontology release archive on AWS S3.
#     The 2026 releases are listed in the bucket but not publicly readable (HTTP 403, checked 2026-11-04),
#     and release.geneontology.org was not reachable from the cloud environment, so 2025-10-10 is the newest release used.
#   HGNC complete set, monthly archive of 2025-10-07 (Ensembl gene -> UniProt accessions, for the airway genes).
#   The airway data package (Himes et al. 2014) from the Bioconductor source mirror on GitHub, at a fixed commit.
# Usage: bash download.sh [data_dir]      (default: ./data); checks md5sums.txt at the end
set -euo pipefail
cd "$(dirname "$0")"
D=${1:-data}
mkdir -p "$D"
GO=https://go-data-product-release.s3.amazonaws.com
for r in 2016-01-01 2020-01-01 2023-01-01 2025-10-10; do
  mkdir -p "$D/$r"
  [ -s "$D/$r/go-basic.obo" ] || curl -fsSL -o "$D/$r/go-basic.obo" "$GO/$r/ontology/go-basic.obo"
  [ -s "$D/$r/goa_human.gaf.gz" ] || curl -fsSL -o "$D/$r/goa_human.gaf.gz" "$GO/$r/annotations/goa_human.gaf.gz"
done
[ -s "$D/hgnc_complete_set_2025-10-07.txt" ] || curl -fsSL -o "$D/hgnc_complete_set_2025-10-07.txt" \
  https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2025-10-07.txt
if [ ! -s "$D/airway.RData" ]; then
  tmp=$(mktemp -d)
  git -C "$tmp" init -q airway
  git -C "$tmp/airway" fetch -q --depth 1 https://github.com/bioconductor-source/airway 596678815f4ade04a9711997150949b9782aeddc
  git -C "$tmp/airway" checkout -q FETCH_HEAD
  cp "$tmp/airway/data/airway.RData" "$D/airway.RData"
  rm -rf "$tmp"
fi
(cd "$D" && md5sum -c "$OLDPWD/md5sums.txt")
