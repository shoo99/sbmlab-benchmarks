#!/usr/bin/env bash
# Downloads what run.sh needs into <data dir>:
#   fastq/<sample>_<rep>_R1/_R2.fastq.gz  the runs listed in runs.tsv, from ENA (lanes concatenated, optionally subsampled)
#   ref/genome.fa, annotation.gtf, transcripts.fa   GENCODE human (GRCh38 primary assembly, comprehensive CHR annotation)
#   ref/ercc.fa, ercc.gtf, ercc_truth.tsv            ERCC92 sequences and the Mix 1 / Mix 2 design (Thermo Fisher)
#
# Usage: ./download.sh <data dir> [read pairs per replicate, 0 = all]
# Environment: GENCODE (release, default 46); ERCC_ZIP_URL and ERCC_CONC_URL (Thermo Fisher links; check they still
# resolve); subsampling needs seqtk and uses seed 11 for both mates, so pairs stay matched.
# Re-running skips files that are already complete.
set -euo pipefail
cd "$(dirname "$0")"
DATA=$(mkdir -p "${1:?usage: ./download.sh <data dir> [read pairs per replicate]}" && cd "$1" && pwd)
N=${2:-0}
GENCODE=${GENCODE:-46}
ERCC_ZIP_URL=${ERCC_ZIP_URL:-https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip}
ERCC_CONC_URL=${ERCC_CONC_URL:-https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt}
mkdir -p "$DATA/ref" "$DATA/fastq" "$DATA/download/runs"
grep -v '^#' runs.tsv | awk -F'\t' '$1 != "sample" && NF >= 3' > "$DATA/download/runs.selected.tsv" || true   # skip header lines
[ -s "$DATA/download/runs.selected.tsv" ] || { echo "runs.tsv lists no runs; fill it in first (see ./list_runs.sh)"; exit 1; }

get() {   # get <url> <file>: download once, atomically
  if [ ! -s "$2" ]; then curl -fL --retry 3 -o "$2.part" "$1" && mv "$2.part" "$2"; fi
}

# ---- references
G=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$GENCODE
get "$G/GRCh38.primary_assembly.genome.fa.gz" "$DATA/download/genome.fa.gz"
get "$G/gencode.v$GENCODE.annotation.gtf.gz" "$DATA/download/annotation.gtf.gz"
get "$G/gencode.v$GENCODE.transcripts.fa.gz" "$DATA/download/transcripts.fa.gz"
for f in genome.fa annotation.gtf transcripts.fa; do
  [ -s "$DATA/ref/$f" ] || zcat "$DATA/download/$f.gz" > "$DATA/ref/$f"
done
echo "GENCODE release $GENCODE" > "$DATA/ref/SOURCE.txt"

get "$ERCC_ZIP_URL" "$DATA/download/ERCC92.zip"
unzip -o -q "$DATA/download/ERCC92.zip" -d "$DATA/download/ercc"
cp "$DATA/download/ercc/ERCC92.fa" "$DATA/ref/ercc.fa"
cp "$DATA/download/ercc/ERCC92.gtf" "$DATA/ref/ercc.gtf"
get "$ERCC_CONC_URL" "$DATA/download/ercc_concentrations.txt"
python3 "$(pwd)/ercc_truth.py" "$DATA/download/ercc_concentrations.txt" "$DATA/ref/ercc_truth.tsv"
echo "ERCC: $ERCC_ZIP_URL, $ERCC_CONC_URL" >> "$DATA/ref/SOURCE.txt"

# ---- reads
while IFS=$'\t' read -r sample rep acc; do
  [[ $sample =~ ^[ABCD]$ ]] || { echo "runs.tsv: sample must be A, B, C or D (got '$sample')"; exit 1; }
  urls=$(curl -fsS --retry 3 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$acc&result=read_run&fields=fastq_ftp&format=tsv" | tail -1 | cut -f2)
  for u in ${urls//;/ }; do
    f=$(basename "$u")
    case $f in *_1.fastq.gz|*_2.fastq.gz|"$acc".fastq.gz) get "https://$u" "$DATA/download/runs/$f" ;; esac
  done
done < "$DATA/download/runs.selected.tsv"

cut -f1,2 "$DATA/download/runs.selected.tsv" | sort -u | while IFS=$'\t' read -r sample rep; do
  accs=$(awk -F'\t' -v s="$sample" -v r="$rep" '$1 == s && $2 == r { print $3 }' "$DATA/download/runs.selected.tsv")
  for mate in 1 2; do
    parts=(); for acc in $accs; do
      if [ -f "$DATA/download/runs/${acc}_$mate.fastq.gz" ]; then parts+=("$DATA/download/runs/${acc}_$mate.fastq.gz")
      elif [ "$mate" = 1 ] && [ -f "$DATA/download/runs/$acc.fastq.gz" ]; then parts+=("$DATA/download/runs/$acc.fastq.gz"); fi
    done
    [ ${#parts[@]} -gt 0 ] || continue
    out="$DATA/fastq/${sample}_${rep}_R$mate.fastq.gz"
    [ -s "$out" ] && continue
    if [ "$N" -gt 0 ]; then
      cat "${parts[@]}" | seqtk sample -s 11 - "$N" | gzip -1 > "$out.part"
    else
      cat "${parts[@]}" > "$out.part"
    fi
    mv "$out.part" "$out"
  done
done
echo "done: $(ls "$DATA"/fastq/*_R1.fastq.gz | wc -l) replicates in $DATA/fastq; references in $DATA/ref"
