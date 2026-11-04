#!/usr/bin/env bash
# Downloads the inputs of the HG002 chr20 variant-calling benchmark into data/ (about 2.2 GB).
#   HG002.chr20.bam(.bai)   HG002 (NA24385) Illumina NovaSeq PCR-free WGS, 35x, BWA-MEM to GRCh38 no-alt, duplicates
#                           marked, chr20 only: the DeepVariant case-study input (Google Cloud Storage, deepvariant bucket)
#   GRCh38_no_alt.fasta.gz  GRCh38 no-alt analysis set (GCA_000001405.15), the reference of that BAM (GIAB S3 mirror)
#   truth.vcf.gz(.tbi), truth.bed   GIAB HG002 v4.2.1 benchmark for GRCh38 (small variants, chr1-22) and its
#                           benchmark regions without inconsistent sites
#   strata/*.bed.gz         GIAB genome stratifications v3.3 for GRCh38: all difficult regions, not in all difficult
#                           regions, tandem repeats + homopolymers (slop 5), segmental duplications, low mappability
# Usage: bash download.sh    (files that already exist are kept; md5 sums go to data/md5.txt)
set -euo pipefail
cd "$(dirname "$0")" || exit 1
mkdir -p data/strata
DV=https://storage.googleapis.com/deepvariant/case-study-testdata
GIAB=https://giab.s3.amazonaws.com/release
TRUTH=$GIAB/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38
STRAT=$GIAB/genome-stratifications/v3.3/GRCh38@all
fetch() {   # fetch <url> <file>
  if [ -s "data/$2" ]; then echo "have $2"; return; fi
  echo "get  $2"
  curl -fsSL --retry 3 -o "data/$2.part" "$1" && mv "data/$2.part" "data/$2"
}
fetch "$DV/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam" HG002.chr20.bam
fetch "$DV/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai" HG002.chr20.bam.bai
fetch "$GIAB/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz" GRCh38_no_alt.fasta.gz
fetch "$TRUTH/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" truth.vcf.gz
fetch "$TRUTH/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" truth.vcf.gz.tbi
fetch "$TRUTH/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" truth.bed
for f in Union/GRCh38_alldifficultregions.bed.gz Union/GRCh38_notinalldifficultregions.bed.gz \
         LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz SegmentalDuplications/GRCh38_segdups.bed.gz \
         Mappability/GRCh38_lowmappabilityall.bed.gz; do
  fetch "$STRAT/$f" "strata/$(basename "$f")"
done
(cd data && md5sum HG002.chr20.bam GRCh38_no_alt.fasta.gz truth.vcf.gz truth.bed strata/*.bed.gz > md5.txt && cat md5.txt)
echo "next: bash run.sh   (see README for the three conda environments)"
