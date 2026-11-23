#!/usr/bin/env bash
# DIA-NN 2.3.2 on the Orbitrap (Q Exactive HF-X) DIA runs of PXD028735: the same six samples as the Orbitrap DDA runs
# searched with MaxQuant and FragPipe (LFQ_Orbitrap_AIF_Condition_{A,B}_Sample_{Alpha,Beta,Gamma}_01, as mzML from
# convert_orbitrap.sh). Despite "AIF" in the names these are DIA runs with 8 m/z isolation windows over 400-1000 m/z,
# staggered by 4 m/z between cycles (151 distinct window centres); DIA-NN handles the overlapping windows itself.
# This removes the instrument difference of run_diann.sh (TripleTOF SWATH) from the DIA-NN vs MaxQuant / FragPipe
# comparison. Same DIA-NN settings as run_diann.sh and the same three predicted libraries; only the raw data change:
#   target   TARGET_LIB  + TARGET_FASTA   (target only: the diann-hye library and hye_sp.fasta; for the targets lost)
#   foreign  FOREIGN_LIB + FOREIGN_FASTA  (target + Arabidopsis entrapment: run_diann.sh's library, db/ from make_db.sh)
#   shuffle  SHUFFLE_LIB + SHUFFLE_FASTA  (target + shuffled entrapment: the same from make_db.sh shuffle)
# Per database: MBR on (--reanalyse; the first-pass .quant files go to <out>/quant_<name>), then MBR off (--use-quant
# reuses those .quant files, no second pass), as for the SWATH runs.
# Usage: DIANN=/path/to/diann TARGET_LIB=... TARGET_FASTA=... FOREIGN_LIB=... FOREIGN_FASTA=... SHUFFLE_LIB=...
#        SHUFFLE_FASTA=... bash run_diann_orbitrap.sh <mzML folder> <out folder> [threads]
# Writes <out>/{target,foreign,shuffle}_mbr_{on,off}/ (report.parquet, diann.log). On the workstation (24 threads,
# 2026-11-18) DIA-NN was the Windows build called from WSL; with a diann.exe the paths are converted with wslpath -w.
set -u
MZML=${1:?mzML folder}; OUT=${2:?out folder}; THREADS=${3:-24}
DIANN=${DIANN:-diann}
win() { if [[ $DIANN == *.exe ]]; then wslpath -w "$1"; else echo "$1"; fi; }
log() { echo "== $(date '+%F %T') $*"; }
mkdir -p "$OUT" || exit 1
F=()
for r in "$MZML"/LFQ_Orbitrap_AIF_Condition_*_01.mzML; do [ -e "$r" ] && F+=(--f "$(win "$r")"); done
[ ${#F[@]} -eq 12 ] || { log "expected 6 mzML files in $MZML, found $(( ${#F[@]} / 2 ))"; exit 1; }
check() {   # check <log> <exit>: DIA-NN exits 0 even when it reads nothing, and a corrupt file kills it with 9
  if [ "$2" != 0 ] || ! grep -q "^6 files will be processed" "$1" || grep -q "^ERROR\|raw file corruption" "$1" \
     || [ ! -s "$(dirname "$1")/report.parquet" ]; then
    log "FAILED: $1 (exit $2)"; grep "^ERROR\|WARNING: possible\|files will be processed" "$1"; exit 1
  fi
}
COMMON=(--threads "$THREADS" --verbose 1 --qvalue 0.01 --matrices --cut 'K*,R*' --missed-cleavages 1 --min-pep-len 7
        --max-pep-len 30 --min-pr-charge 2 --max-pr-charge 4 --min-pr-mz 300 --max-pr-mz 1800 --unimod4 --var-mods 1
        --var-mod UniMod:35,15.994915,M --rt-profiling)

run_db() {   # run_db <name> <lib> <fasta>
  local n=$1 lib=$2 fa=$3 rc
  mkdir -p "$OUT/quant_$n" "$OUT/${n}_mbr_on" "$OUT/${n}_mbr_off"
  log "$n MBR on"
  "$DIANN" "${F[@]}" --lib "$(win "$lib")" --fasta "$(win "$fa")" --temp "$(win "$OUT/quant_$n")" \
    --out "$(win "$OUT/${n}_mbr_on/report.parquet")" "${COMMON[@]}" --reanalyse > "$OUT/${n}_mbr_on/diann.log" 2>&1
  rc=$?; log "$n MBR on exit $rc"; check "$OUT/${n}_mbr_on/diann.log" $rc
  "$DIANN" "${F[@]}" --lib "$(win "$lib")" --fasta "$(win "$fa")" --temp "$(win "$OUT/quant_$n")" \
    --out "$(win "$OUT/${n}_mbr_off/report.parquet")" "${COMMON[@]}" --use-quant > "$OUT/${n}_mbr_off/diann.log" 2>&1
  rc=$?; log "$n MBR off exit $rc"; check "$OUT/${n}_mbr_off/diann.log" $rc
}

run_db target  "${TARGET_LIB:?}"  "${TARGET_FASTA:?}"
run_db foreign "${FOREIGN_LIB:?}" "${FOREIGN_FASTA:?}"
run_db shuffle "${SHUFFLE_LIB:?}" "${SHUFFLE_FASTA:?}"
log "done"
