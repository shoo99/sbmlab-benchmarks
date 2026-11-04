#!/usr/bin/env bash
# HG002 chr20 germline small-variant calling: bcftools, GATK HaplotypeCaller and DeepVariant on the same BAM, scored
# with hap.py (vcfeval engine) against the GIAB v4.2.1 benchmark, overall and per GIAB stratification (issue #14).
#
# Usage: bash run.sh [step ...]     steps: prepare bcftools gatk deepvariant happy   (default: all, in that order)
#   env: VC_BIN     bin folder with samtools, bcftools, gatk, rtg      (conda env from envs/vc.yml)
#        DV_RUN     command that runs DeepVariant's run_deepvariant (default: the official Docker image,
#                   docker run google/deepvariant:1.10.0 with data/ and work/ mounted at the same paths)
#        HAPPY_BIN  bin folder with hap.py and rtg                      (conda env from envs/happy.yml)
#        THREADS    (default 4)   DATA (default data)   WORK (default work)   COVERAGE (default 35; e.g. 15 downsamples)
# Writes: $WORK/<caller>.vcf.gz (not tracked), happy/<caller>_<coverage>x.extended.csv (trimmed, tracked),
#         timings.tsv (tracked: caller, coverage, wall seconds, CPU seconds, peak RSS of the largest process in MB).
# Then: python3 score.py > score.out && python3 figure.py
set -euo pipefail
cd "$(dirname "$0")" || exit 1
THREADS=${THREADS:-4}; DATA=${DATA:-data}; WORK=${WORK:-work}; COV=${COVERAGE:-35}
VC_BIN=${VC_BIN:-}; HAPPY_BIN=${HAPPY_BIN:-}
tool() { if [ -n "$2" ]; then echo "$2/$1"; else echo "$1"; fi; }
SAMTOOLS=$(tool samtools "$VC_BIN"); BCFTOOLS=$(tool bcftools "$VC_BIN"); GATK=$(tool gatk "$VC_BIN"); RTG=$(tool rtg "$VC_BIN")
HAPPY=$(tool hap.py "$HAPPY_BIN")
REF=$DATA/GRCh38_no_alt.fasta; REGION=chr20
STEPS=("$@"); [ ${#STEPS[@]} -eq 0 ] && STEPS=(prepare bcftools gatk deepvariant happy)
mkdir -p "$WORK" happy
[ -s timings.tsv ] || printf 'caller\tcoverage\twall_s\tcpu_s\tpeak_rss_mb\n' > timings.tsv

timed() {   # timed <caller label> <command ...>: run, append wall / CPU time and peak RSS to timings.tsv
  local label=$1; shift
  local log="$WORK/$label.time"
  /usr/bin/time -v -o "$log" "$@"
  local wall cpu rss
  wall=$(awk -F': ' '/Elapsed \(wall clock\)/ {n = split($2, a, ":"); s = 0; for (i = 1; i <= n; i++) s = s * 60 + a[i]; print s}' "$log")
  cpu=$(awk -F': ' '/User time/ {u = $2} /System time/ {s = $2} END {print u + s}' "$log")
  rss=$(awk -F': ' '/Maximum resident set size/ {printf "%.0f", $2 / 1024}' "$log")
  grep -v "^$label	$COV	" timings.tsv > timings.tmp || true; mv timings.tmp timings.tsv
  printf '%s\t%s\t%.0f\t%.0f\t%s\n' "$label" "$COV" "$wall" "$cpu" "$rss" >> timings.tsv
}

bam() {     # the input BAM, downsampled from 35x when COVERAGE < 35 (seed 7)
  if [ "$COV" = 35 ]; then echo "$DATA/HG002.chr20.bam"; return; fi
  local out="$WORK/HG002.chr20.${COV}x.bam"
  if [ ! -s "$out" ]; then
    frac=$(awk -v c="$COV" 'BEGIN {printf "%.4f", c / 35}')
    "$SAMTOOLS" view -@ "$THREADS" -b -s "7${frac#0}" -o "$out" "$DATA/HG002.chr20.bam" && "$SAMTOOLS" index "$out"
  fi
  echo "$out"
}

for step in "${STEPS[@]}"; do
  case $step in
  prepare)
    # reference: plain FASTA + .fai + .dict; truth and benchmark regions for chr20; stratification table for hap.py
    [ -s "$REF" ] || gzip -dc "$DATA/GRCh38_no_alt.fasta.gz" > "$REF"
    [ -s "$REF.fai" ] || "$SAMTOOLS" faidx "$REF"
    [ -s "${REF%.fasta}.dict" ] || "$GATK" CreateSequenceDictionary -R "$REF" -O "${REF%.fasta}.dict" > "$WORK/dict.log" 2>&1
    [ -d "$DATA/chr20.sdf" ] || { "$SAMTOOLS" faidx "$REF" chr20 > "$WORK/chr20.fa" && "$RTG" format -o "$DATA/chr20.sdf" "$WORK/chr20.fa" > /dev/null; }
    "$BCFTOOLS" view -r chr20 -Oz -o "$DATA/truth.chr20.vcf.gz" "$DATA/truth.vcf.gz" && "$BCFTOOLS" index -t -f "$DATA/truth.chr20.vcf.gz"
    awk '$1 == "chr20"' "$DATA/truth.bed" > "$DATA/truth.chr20.bed"
    mkdir -p "$DATA/strata.chr20"
    : > "$DATA/strata.tsv"
    for f in "$DATA"/strata/*.bed.gz; do
      name=$(basename "$f" .bed.gz | sed 's/^GRCh38_//')
      gzip -dc "$f" | awk '$1 == "chr20"' > "$DATA/strata.chr20/$name.bed"
      printf '%s\t%s\n' "$name" "$(cd "$DATA/strata.chr20" && pwd)/$name.bed" >> "$DATA/strata.tsv"
    done
    echo "prepared: $(grep -vc '^#' <(gzip -dc "$DATA/truth.chr20.vcf.gz")) truth variants on chr20, $(wc -l < "$DATA/truth.chr20.bed") benchmark intervals, $(wc -l < "$DATA/strata.tsv") strata"
    ;;
  bcftools)
    # mpileup (MAPQ >= 20, base quality >= 20, depth cap 250) | call -m (multiallelic caller), variants only; a soft
    # filter LowQual (QUAL < 20) keeps every call so hap.py scores both all calls and the filtered set
    B=$(bam); out="$WORK/bcftools_${COV}x.vcf.gz"
    timed bcftools bash -c "'$BCFTOOLS' mpileup --threads $THREADS -f '$REF' -r $REGION -q 20 -Q 20 -d 250 -a AD,DP -Ou '$B' \
      | '$BCFTOOLS' call --threads $THREADS -mv -Ou | '$BCFTOOLS' filter -s LowQual -e 'QUAL<20' -Oz -o '$out'"
    "$BCFTOOLS" index -t -f "$out"
    ;;
  gatk)
    # HaplotypeCaller, chr20 split into THREADS equal pieces run in parallel, joined with bcftools concat; GATK's
    # recommended hard filters as FILTER flags (SNP: QD < 2, FS > 60, MQ < 40, MQRankSum < -12.5, ReadPosRankSum < -8,
    # SOR > 3; indel: QD < 2, FS > 200, ReadPosRankSum < -20, SOR > 10). VQSR is not run: it needs many samples or a
    # whole genome of variants and external resource panels.
    B=$(bam); out="$WORK/gatk_${COV}x.vcf.gz"
    len=$(awk '$1 == "chr20" {print $2}' "$REF.fai"); step_bp=$(( (len + THREADS - 1) / THREADS ))
    pieces=(); for i in $(seq 0 $((THREADS - 1))); do s=$((i * step_bp + 1)); e=$(( (i + 1) * step_bp )); [ $e -gt $len ] && e=$len; pieces+=("chr20:$s-$e"); done
    cat > "$WORK/gatk_scatter.sh" <<SH
set -e
i=0
for iv in ${pieces[*]}; do
  '$GATK' --java-options "-Xmx3g -XX:ParallelGCThreads=1" HaplotypeCaller -R '$REF' -I '$B' -L \$iv -O '$WORK'/gatk_part\$i.vcf.gz > '$WORK'/gatk_part\$i.log 2>&1 &
  i=\$((i + 1))
done
wait
SH
    timed gatk bash "$WORK/gatk_scatter.sh"
    parts=(); for i in $(seq 0 $((THREADS - 1))); do parts+=("$WORK/gatk_part$i.vcf.gz"); done
    "$BCFTOOLS" concat -a -D -Oz -o "$WORK/gatk_raw_${COV}x.vcf.gz" "${parts[@]}" 2> /dev/null || "$BCFTOOLS" concat -Oz -o "$WORK/gatk_raw_${COV}x.vcf.gz" "${parts[@]}"
    "$BCFTOOLS" index -t -f "$WORK/gatk_raw_${COV}x.vcf.gz"
    "$GATK" SelectVariants -V "$WORK/gatk_raw_${COV}x.vcf.gz" -select-type SNP -O "$WORK/gatk_snp.vcf.gz" > /dev/null 2>&1
    "$GATK" SelectVariants -V "$WORK/gatk_raw_${COV}x.vcf.gz" -select-type INDEL -select-type MIXED -O "$WORK/gatk_indel.vcf.gz" > /dev/null 2>&1
    "$GATK" VariantFiltration -V "$WORK/gatk_snp.vcf.gz" -O "$WORK/gatk_snp.f.vcf.gz" \
      -filter "QD < 2.0" --filter-name QD2 -filter "FS > 60.0" --filter-name FS60 -filter "MQ < 40.0" --filter-name MQ40 \
      -filter "MQRankSum < -12.5" --filter-name MQRankSum-12.5 -filter "ReadPosRankSum < -8.0" --filter-name ReadPosRankSum-8 \
      -filter "SOR > 3.0" --filter-name SOR3 > /dev/null 2>&1
    "$GATK" VariantFiltration -V "$WORK/gatk_indel.vcf.gz" -O "$WORK/gatk_indel.f.vcf.gz" \
      -filter "QD < 2.0" --filter-name QD2 -filter "FS > 200.0" --filter-name FS200 \
      -filter "ReadPosRankSum < -20.0" --filter-name ReadPosRankSum-20 -filter "SOR > 10.0" --filter-name SOR10 > /dev/null 2>&1
    "$BCFTOOLS" concat -a -Oz -o "$out" "$WORK/gatk_snp.f.vcf.gz" "$WORK/gatk_indel.f.vcf.gz" && "$BCFTOOLS" index -t -f "$out"
    ;;
  deepvariant)
    # DeepVariant 1.10.0, WGS model (with its small model for easy sites), THREADS shards, via the official Docker
    # image unless DV_RUN names another way to call run_deepvariant (e.g. singularity exec ... run_deepvariant).
    # RefCall records carry FILTER RefCall, so hap.py's PASS rows are its calls.
    B=$(bam); out="$WORK/deepvariant_${COV}x.vcf.gz"
    D=$(cd "$DATA" && pwd -P); W=$(cd "$WORK" && pwd -P)
    DV_RUN=${DV_RUN:-docker run --rm -v "$D:$D" -v "$W:$W" -w "$W" google/deepvariant:1.10.0 /opt/deepvariant/bin/run_deepvariant}
    timed deepvariant bash -c "$DV_RUN --model_type=WGS --ref='$D/$(basename "$REF")' --reads='$(cd "$(dirname "$B")" && pwd -P)/$(basename "$B")' \
      --regions=$REGION --output_vcf='$W/$(basename "$out")' --num_shards=$THREADS --intermediate_results_dir='$W/dv_tmp_${COV}x'"
    ;;
  happy)
    # hap.py with the vcfeval engine against the chr20 truth, in the benchmark regions, with the GIAB strata
    for caller in bcftools gatk deepvariant; do
      q="$WORK/${caller}_${COV}x.vcf.gz"; [ -s "$q" ] || { echo "skip $caller (no $q)"; continue; }
      p="$WORK/happy_${caller}_${COV}x"
      HGREF="$REF" "$HAPPY" "$DATA/truth.chr20.vcf.gz" "$q" -f "$DATA/truth.chr20.bed" -r "$REF" -l chr20 \
        --engine vcfeval --engine-vcfeval-template "$DATA/chr20.sdf" --stratification "$DATA/strata.tsv" \
        --threads "$THREADS" -o "$p" > "$p.log" 2>&1
      # keep the rows scored here: Subtype *, Filter ALL / PASS, SNP / INDEL, overall and per stratum
      python3 -c '
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
keep = ["Type", "Subset", "Filter", "METRIC.Recall", "METRIC.Precision", "METRIC.F1_Score", "TRUTH.TOTAL", "TRUTH.TP",
        "TRUTH.FN", "QUERY.TOTAL", "QUERY.FP", "QUERY.UNK", "FP.gt", "FP.al", "Subset.Size", "Subset.IS_CONF.Size"]
w = csv.DictWriter(open(sys.argv[2], "w", newline=""), fieldnames=keep, extrasaction="ignore", lineterminator="\n")
w.writeheader()
for r in rows:
    if r["Subtype"] == "*" and r["Filter"] in ("ALL", "PASS") and r["Type"] in ("SNP", "INDEL"):
        w.writerow(r)
' "$p.extended.csv" "happy/${caller}_${COV}x.extended.csv"
      echo "hap.py $caller ${COV}x: $(grep -c . "happy/${caller}_${COV}x.extended.csv") rows"
    done
    ;;
  *) echo "unknown step $step"; exit 1 ;;
  esac
done
