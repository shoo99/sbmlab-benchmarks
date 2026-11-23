#!/usr/bin/env bash
# .raw -> .mzML for the DIA-NN Orbitrap DIA searches (run_diann_orbitrap.sh). DIA-NN 2.3.2 reads Thermo .raw only with
# .NET 8, which was not installed on the workstation, so the six runs were converted first.
# As run (2026-11-18): ProteoWizard msconvert 3.0.25062 (Windows build, called from WSL), --mzML --32 and vendor peak
# picking as the first filter (--filter "peakPicking vendor msLevel=1-"), following the DIA-NN README (32-bit, index,
# vendor peak picking first). The README also lists no compression, but this msconvert version writes zlib
# compression unless told otherwise, and the command had no switch for it: the six mzML files used carry
# "zlib compression" and "32-bit float" on every binary array. zlib is lossless, so the spectra are the same.
# Usage: MSCONVERT=/path/to/msconvert bash convert_orbitrap.sh <folder with the .raw files> <mzML folder> <file stems...>
#   e.g. stems LFQ_Orbitrap_AIF_Condition_{A,B}_Sample_{Alpha,Beta,Gamma}_01 (3 conversions in parallel).
# With a Windows msconvert.exe under WSL, the paths given to it are converted with wslpath -w.
set -u
RAW=${1:?folder with the .raw files}; OUT=${2:?mzML folder}; shift 2
MSCONVERT=${MSCONVERT:-msconvert}
mkdir -p "$OUT"
win() { if [[ $MSCONVERT == *.exe ]]; then wslpath -w "$1"; else echo "$1"; fi; }
conv() {   # skips a file whose mzML is complete (ends with the index)
  local stem=$1
  if [ -s "$OUT/$stem.mzML" ] && tail -c 200 "$OUT/$stem.mzML" | grep -q "</indexedmzML>"; then echo "HAVE $stem"; return; fi
  "$MSCONVERT" "$(win "$RAW/$stem.raw")" --mzML --32 --filter "peakPicking vendor msLevel=1-" -o "$(win "$OUT")" \
    > "$OUT/$stem.conv.log" 2>&1
  tail -c 200 "$OUT/$stem.mzML" | grep -q "</indexedmzML>" && echo "OK $stem" || echo "FAIL $stem"
}
export -f conv win; export RAW OUT MSCONVERT
printf '%s\n' "$@" | xargs -P 3 -I{} bash -c 'conv {}'
echo CONVDONE
