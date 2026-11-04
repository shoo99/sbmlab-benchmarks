#!/usr/bin/env bash
# Taxonomic profilers on the ZymoBIOMICS mock communities (issue #15). Profiles go to results/, scored by score.py.
#
# Usage: bash run.sh <step> [dataset ...]        datasets: D6300 D6310 (default both)
#   kraken   Kraken2 + Bracken (species level) with every database under db/ (k2_standard_08gb, k2_pluspf_08gb, ...)
#   metaphlan, motus, sylph   the profilers whose databases are too large or not reachable from the cloud environment;
#            run on a workstation with their own conda environments (see README) -> results/<dataset>_<tool>.profile.tsv
#   env: THREADS (4), N_SPOTS (all; e.g. 2000000 keeps the first 2 M read pairs), MG_BIN (bin folder with kraken2,
#        bracken, fastq-dump), METAPHLAN_DB, MOTUS_DB, SYLPH_DB
# Reads: both mates of every pair are classified as single reads, streamed from the .sra with fastq-dump (no FASTQ on
# disk). Bracken: species level (-l S), at least 10 reads (-t 10), read length 150 (D6300) or 100 (D6310).
# Timing: wall and CPU time and peak RSS of the largest process go to timings.tsv.
set -euo pipefail
cd "$(dirname "$0")" || exit 1
STEP=${1:?step: kraken | metaphlan | motus | sylph}; shift
SETS=("$@"); [ ${#SETS[@]} -eq 0 ] && SETS=(D6300 D6310)
THREADS=${THREADS:-4}; N_SPOTS=${N_SPOTS:-}
B=${MG_BIN:+$MG_BIN/}
declare -A RUN=([D6300]=ERR2984773 [D6310]=ERR2935805) RLEN=([D6300]=150 [D6310]=100)
mkdir -p results work
[ -s timings.tsv ] || printf 'dataset\ttool\twall_s\tcpu_s\tpeak_rss_mb\n' > timings.tsv

timed() {   # timed <dataset> <tool label> <command ...>
  local ds=$1 label=$2; shift 2
  local log="work/${ds}_${label}.time"
  /usr/bin/time -v -o "$log" "$@"
  local wall cpu rss
  wall=$(awk -F': ' '/Elapsed \(wall clock\)/ {n = split($2, a, ":"); s = 0; for (i = 1; i <= n; i++) s = s * 60 + a[i]; print s}' "$log")
  cpu=$(awk -F': ' '/User time/ {u = $2} /System time/ {s = $2} END {print u + s}' "$log")
  rss=$(awk -F': ' '/Maximum resident set size/ {printf "%.0f", $2 / 1024}' "$log")
  grep -v "^$ds	$label	" timings.tsv > timings.tmp || true; mv timings.tmp timings.tsv
  printf '%s\t%s\t%.0f\t%.0f\t%s\n' "$ds" "$label" "$wall" "$cpu" "$rss" >> timings.tsv
}
reads() {   # reads <dataset>: FASTQ of both mates on stdout
  "${B}fastq-dump" --stdout --split-spot ${N_SPOTS:+-X $N_SPOTS} "data/${RUN[$1]}.sra" 2> /dev/null
}

for ds in "${SETS[@]}"; do
  case $STEP in
  kraken)
    for db in db/*/; do
      db=${db%/}; name=$(basename "$db"); [ -s "$db/hash.k2d" ] || continue
      rep="results/${ds}_${name}.kreport"
      timed "$ds" "kraken2 ($name)" bash -c "$(declare -f reads); B='$B'; N_SPOTS='$N_SPOTS'; declare -A RUN=([D6300]=ERR2984773 [D6310]=ERR2935805);
        reads $ds | '${B}kraken2' --db '$db' --threads $THREADS --report '$rep' --output /dev/null /dev/stdin 2> 'work/${ds}_${name}.kraken.log'"
      "${B}bracken" -d "$db" -i "$rep" -o "results/${ds}_${name}.bracken" -w "work/${ds}_${name}.breport" \
        -r "${RLEN[$ds]}" -l S -t 10 > "work/${ds}_${name}.bracken.log" 2>&1
      echo "$ds $name: $(grep -E 'processed|classified' "work/${ds}_${name}.kraken.log" | tr '\n' ' ')"
    done
    ;;
  metaphlan)
    # MetaPhlAn 4 (bowtie2 database, about 20 GB): relative abundances of species (s__) rows
    timed "$ds" metaphlan bash -c "$(declare -f reads); B='$B'; N_SPOTS='$N_SPOTS'; declare -A RUN=([D6300]=ERR2984773 [D6310]=ERR2935805);
      reads $ds > work/$ds.fastq && metaphlan work/$ds.fastq --input_type fastq --nproc $THREADS ${METAPHLAN_DB:+--db_dir $METAPHLAN_DB} \
      --bowtie2out work/$ds.bt2.bz2 -o work/${ds}_metaphlan.txt && rm -f work/$ds.fastq"
    awk -F'\t' 'BEGIN {OFS = "\t"; print "taxid", "species", "percent"}
      $1 ~ /\|s__[^|]*$/ {n = split($1, a, "s__"); sp = a[n]; gsub("_", " ", sp); split($2, t, "|"); print t[length(t)], sp, $3}' \
      "work/${ds}_metaphlan.txt" > "results/${ds}_metaphlan.profile.tsv"
    ;;
  motus)
    # mOTUs 3 (database about 3.5 GB): species-level relative abundances (-k mOTU is too fine; -k species)
    timed "$ds" motus bash -c "$(declare -f reads); B='$B'; N_SPOTS='$N_SPOTS'; declare -A RUN=([D6300]=ERR2984773 [D6310]=ERR2935805);
      reads $ds > work/$ds.fastq && motus profile -s work/$ds.fastq -t $THREADS -k species -c ${MOTUS_DB:+-db $MOTUS_DB} \
      -o work/${ds}_motus.txt && rm -f work/$ds.fastq"
    awk -F'\t' 'BEGIN {OFS = "\t"; print "taxid", "species", "percent"; s = 0}
      !/^#/ && $1 != "unassigned" {n[$1] = $2; s += $2} END {for (k in n) if (n[k] > 0) print "", k, 100 * n[k] / s}' \
      "work/${ds}_motus.txt" > "results/${ds}_motus.profile.tsv"
    ;;
  sylph)
    # sylph with the GTDB r220 database (about 14 GB) and sylph-tax for species names
    timed "$ds" sylph bash -c "$(declare -f reads); B='$B'; N_SPOTS='$N_SPOTS'; declare -A RUN=([D6300]=ERR2984773 [D6310]=ERR2935805);
      reads $ds > work/$ds.fastq && sylph profile '${SYLPH_DB:?set SYLPH_DB}' work/$ds.fastq -t $THREADS -o work/${ds}_sylph.tsv \
      && rm -f work/$ds.fastq"
    sylph-tax taxprof "work/${ds}_sylph.tsv" -t GTDB_r220 -o "work/${ds}_sylph_"
    awk -F'\t' 'BEGIN {OFS = "\t"; print "taxid", "species", "percent"}
      $1 ~ /\|s__[^|]*$/ {n = split($1, a, "s__"); print "", a[n], $2}' "work/${ds}_sylph_"*.sylphmpa > "results/${ds}_sylph.profile.tsv"
    ;;
  *) echo "unknown step $STEP"; exit 1 ;;
  esac
done
