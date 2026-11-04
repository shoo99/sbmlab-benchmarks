#!/usr/bin/env bash
# Runs the whole doublet study on the downloads: prepare.R -> R tools + Scrublet per data set -> collect.py -> score.py
# Usage: bash run.sh [dataset ...]     (default: kang_ctrl kang_stim hgmm_12k; run download.sh first)
#   env: DBL_RSCRIPT, DBL_PYTHON (default Rscript, python3 of the environment built from conda-linux-64.lock)
# Writes work/ (inputs and per-tool scores; not tracked), then cell_scores.tsv.gz, datasets.tsv, scores.csv and the
# printouts prepare.out, tools.out, score.out next to this script. Timings per tool are in tools.out.
set -euo pipefail
cd "$(dirname "$0")" || exit 1
RS=${DBL_RSCRIPT:-Rscript}; PY=${DBL_PYTHON:-python3}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$RS" prepare.R data work "$@" > prepare.out
: > tools.out
for d in work/*/; do
  rc=0; "$RS" run_r_tools.R "$d" >> tools.out 2>&1 || rc=$?
  [ "$rc" -eq 0 ] || { echo "run_r_tools.R failed on $d (exit $rc); see tools.out"; exit 1; }
  rc=0; "$PY" run_scrublet.py "$d" >> tools.out 2>&1 || rc=$?
  [ "$rc" -eq 0 ] || { echo "run_scrublet.py failed on $d (exit $rc); see tools.out"; exit 1; }
done
"$PY" collect.py work .
"$PY" score.py > score.out
cat prepare.out score.out
