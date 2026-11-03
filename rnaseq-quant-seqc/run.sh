#!/usr/bin/env bash
# Quantify every sample with Salmon, kallisto, STAR + featureCounts and RSEM (on STAR's transcriptome alignments).
#
# Usage: ./run.sh <data dir> [threads]
#   <data dir>/ref:   genome.fa, annotation.gtf, transcripts.fa (GENCODE headers are fine), ercc.fa, ercc.gtf
#   <data dir>/fastq: <sample>_<rep>_R1.fastq.gz (+ _R2 for paired-end), sample = A, B, C or D
# Writes <data dir>/work/{salmon,kallisto,star,featurecounts,rsem}/..., samples.txt, versions.tsv, timings.tsv.
# Environment: STRAND = none | forward | reverse (default none; check salmon's lib_format_counts.json for the real
# library type and re-run featureCounts / kallisto / RSEM with it if it is stranded), SALMON_DECOYS=1 for a
# decoy-aware Salmon index (genome as decoys).
set -euo pipefail
cd "$(dirname "$0")"
DATA=$(cd "${1:?usage: ./run.sh <data dir> [threads]}" && pwd)
T=${2:-8}
STRAND=${STRAND:-none}
REF=$DATA/ref
W=$DATA/work
mkdir -p "$W"

case $STRAND in
  none)    FC_S=0; KAL_S=();               RSEM_S=none ;;
  forward) FC_S=1; KAL_S=(--fr-stranded);  RSEM_S=forward ;;
  reverse) FC_S=2; KAL_S=(--rf-stranded);  RSEM_S=reverse ;;
  *) echo "STRAND must be none, forward or reverse"; exit 1 ;;
esac

samples=()
for r1 in "$DATA"/fastq/*_R1.fastq.gz; do
  s=$(basename "$r1" _R1.fastq.gz)
  [[ $s =~ ^[ABCD]_[0-9]+$ ]] || { echo "unexpected fastq name: $r1 (want <A|B|C|D>_<rep>_R1.fastq.gz)"; exit 1; }
  samples+=("$s")
done
[ ${#samples[@]} -gt 0 ] || { echo "no fastq files in $DATA/fastq"; exit 1; }
printf '%s\n' "${samples[@]}" > "$W/samples.txt"
paired() { [ -f "$DATA/fastq/$1_R2.fastq.gz" ]; }
reads() { if paired "$1"; then echo "$DATA/fastq/$1_R1.fastq.gz $DATA/fastq/$1_R2.fastq.gz"; else echo "$DATA/fastq/$1_R1.fastq.gz"; fi; }

echo -e "tool\tversion" > "$W/versions.tsv"
{
  echo -e "salmon\t$(salmon --version 2>&1 | awk '{print $2}')"
  echo -e "kallisto\t$(kallisto version | awk '{print $NF}')"
  echo -e "STAR\t$(STAR --version)"
  echo -e "featureCounts\t$(featureCounts -v 2>&1 | grep -o 'v[0-9][0-9.]*' | head -1 | tr -d v)"
  echo -e "RSEM\t$(rsem-calculate-expression --version | grep -o 'v[0-9][0-9.]*' | head -1 | tr -d v)"
} >> "$W/versions.tsv"
echo -e "tool\tstep\tseconds" > "$W/timings.tsv"
timed() { local tool=$1 step=$2; shift 2; local t0=$SECONDS; "$@"; echo -e "$tool\t$step\t$((SECONDS - t0))" >> "$W/timings.tsv"; }

# shared references: genome + ERCC, annotation + ERCC, transcripts with the ID only (GENCODE headers cut at '|') + ERCC
cat "$REF/genome.fa" "$REF/ercc.fa" > "$W/genome_ercc.fa"
cat "$REF/annotation.gtf" "$REF/ercc.gtf" > "$W/annotation_ercc.gtf"
{ sed 's/|.*//' "$REF/transcripts.fa"; cat "$REF/ercc.fa"; } > "$W/transcripts_ercc.fa"

# ---- Salmon
if [ "${SALMON_DECOYS:-0}" = 1 ]; then
  grep '^>' "$W/genome_ercc.fa" | cut -c2- | cut -d' ' -f1 > "$W/decoys.txt"
  cat "$W/transcripts_ercc.fa" "$W/genome_ercc.fa" > "$W/gentrome.fa"
  timed salmon index salmon index -t "$W/gentrome.fa" -d "$W/decoys.txt" -i "$W/salmon_idx" -p "$T" > "$W/salmon_index.log" 2>&1
else
  timed salmon index salmon index -t "$W/transcripts_ercc.fa" -i "$W/salmon_idx" -p "$T" > "$W/salmon_index.log" 2>&1
fi
for s in "${samples[@]}"; do
  if paired "$s"; then in=(-1 "$DATA/fastq/${s}_R1.fastq.gz" -2 "$DATA/fastq/${s}_R2.fastq.gz"); else in=(-r "$DATA/fastq/${s}_R1.fastq.gz"); fi
  timed salmon quant salmon quant -i "$W/salmon_idx" -l A "${in[@]}" -p "$T" -o "$W/salmon/$s" > "$W/salmon_$s.log" 2>&1
done

# ---- kallisto
timed kallisto index kallisto index -i "$W/kallisto.idx" "$W/transcripts_ercc.fa" > "$W/kallisto_index.log" 2>&1
for s in "${samples[@]}"; do
  if paired "$s"; then se=(); else se=(--single -l 200 -s 30); fi
  mkdir -p "$W/kallisto/$s"                                                  # kallisto does not create parents
  # shellcheck disable=SC2046
  timed kallisto quant kallisto quant -i "$W/kallisto.idx" -o "$W/kallisto/$s" -t "$T" "${KAL_S[@]}" "${se[@]}" $(reads "$s") > "$W/kallisto_$s.log" 2>&1
done

# ---- STAR (genome BAM for featureCounts, transcriptome BAM for RSEM)
glen=$(grep -v '^>' "$W/genome_ercc.fa" | tr -d '\n' | wc -c)
sa=$(awk -v g="$glen" 'BEGIN { n = int(log(g) / log(2) / 2 - 1); print (n < 14 ? n : 14) }')   # STAR's rule for small genomes
mkdir -p "$W/star_idx"
timed STAR index STAR --runMode genomeGenerate --genomeDir "$W/star_idx" --genomeFastaFiles "$W/genome_ercc.fa" \
  --sjdbGTFfile "$W/annotation_ercc.gtf" --genomeSAindexNbases "$sa" --runThreadN "$T" --outFileNamePrefix "$W/star_idx/" > "$W/star_index.log" 2>&1
for s in "${samples[@]}"; do
  mkdir -p "$W/star/$s"
  # shellcheck disable=SC2046
  timed STAR align STAR --genomeDir "$W/star_idx" --readFilesIn $(reads "$s") --readFilesCommand zcat --runThreadN "$T" \
    --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --outFileNamePrefix "$W/star/$s/" > "$W/star_$s.log" 2>&1
done

# ---- featureCounts (gene level, all samples in one table)
mkdir -p "$W/featurecounts"
bams=(); for s in "${samples[@]}"; do bams+=("$W/star/$s/Aligned.out.bam"); done
if paired "${samples[0]}"; then pe=(-p --countReadPairs); else pe=(); fi
timed featureCounts count featureCounts -T "$T" "${pe[@]}" -s "$FC_S" -a "$W/annotation_ercc.gtf" -o "$W/featurecounts/counts.txt" "${bams[@]}" > "$W/featurecounts.log" 2>&1

# ---- RSEM on STAR's transcriptome alignments
mkdir -p "$W/rsem_ref" "$W/rsem"
timed RSEM index rsem-prepare-reference --gtf "$W/annotation_ercc.gtf" "$W/genome_ercc.fa" "$W/rsem_ref/ref" > "$W/rsem_index.log" 2>&1
for s in "${samples[@]}"; do
  if paired "$s"; then pe=(--paired-end); else pe=(); fi
  timed RSEM quant rsem-calculate-expression --alignments "${pe[@]}" --strandedness "$RSEM_S" --no-bam-output -p "$T" \
    "$W/star/$s/Aligned.toTranscriptome.out.bam" "$W/rsem_ref/ref" "$W/rsem/$s" > "$W/rsem_$s.log" 2>&1
done

echo "done: ${#samples[@]} samples; outputs in $W"
