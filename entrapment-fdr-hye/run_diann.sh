#!/usr/bin/env bash
# DIA-NN library-free search of the 6 PXD028735 SWATH runs against the target + entrapment database, with the settings
# of diann-hye (DIA-NN 2.3.2): predicted library, trypsin, 1 missed cleavage, 7-30 aa, charge 2-4, m/z 300-1800,
# fixed carbamidomethyl C, variable Met oxidation (max 1), automatic mass accuracy, 1% FDR, MBR on and off.
# Usage: DIANN=/path/to/diann bash run_diann.sh <folder with the 6 .wiff/.wiff.scan files> [threads]
# Writes searches/diann_mbr_on/ and searches/diann_mbr_off/ (report.parquet, logs). Expect about 6-9 h per search on
# 24 threads (diann-hye: 6 h 17 min with MBR). For the target-only baseline, run it again with FASTA=db/target.fasta and
# OUT_TAG=target_only_. QVALUE (default 0.01, as diann-hye) widens the report for the FDP curve above 1%.
set -euo pipefail
cd "$(dirname "$0")" || exit 1
RAW=${1:?folder with the .wiff files}; THREADS=${2:-24}
DIANN=${DIANN:-diann}; FASTA=${FASTA:-db/target_plus_entrapment.fasta}; TAG=${OUT_TAG:-}; QVALUE=${QVALUE:-0.01}
files=(); for f in "$RAW"/*.wiff; do files+=(--f "$f"); done
[ ${#files[@]} -eq 12 ] || [ ${#files[@]} -eq 6 ] || echo "note: ${#files[@]} / 2 = $(( ${#files[@]} / 2 )) .wiff files found (expected 6)"
lib="searches/${TAG}lib/predicted"
mkdir -p "searches/${TAG}lib"
common=(--threads "$THREADS" --verbose 1 --qvalue "$QVALUE" --matrices --cut 'K*,R*' --missed-cleavages 1 --min-pep-len 7
        --max-pep-len 30 --min-pr-charge 2 --max-pr-charge 4 --min-pr-mz 300 --max-pr-mz 1800 --unimod4 --var-mods 1
        --var-mod UniMod:35,15.994915,M --rt-profiling)
speclib() { ls -t "searches/${TAG}lib/"*.speclib 2>/dev/null | head -1; }
[ -n "$(speclib)" ] || "$DIANN" --fasta "$FASTA" --fasta-search --predictor --gen-spec-lib --out-lib "$lib" "${common[@]}" \
  > "searches/${TAG}lib/libgen.log" 2>&1
LIB=$(speclib); [ -n "$LIB" ] || { echo "no predicted library in searches/${TAG}lib (see libgen.log)"; exit 1; }
for mbr in on off; do
  out="searches/${TAG}diann_mbr_${mbr}"; mkdir -p "$out"
  extra=(); [ "$mbr" = on ] && extra=(--reanalyse)
  rc=0; "$DIANN" "${files[@]}" --lib "$LIB" --fasta "$FASTA" --out "$out/report.parquet" \
    "${common[@]}" "${extra[@]}" > "$out/diann.log" 2>&1 || rc=$?
  echo "$(date '+%F %T') DIA-NN MBR $mbr exit $rc"
done
