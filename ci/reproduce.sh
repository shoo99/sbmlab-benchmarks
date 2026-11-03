#!/usr/bin/env bash
# Re-runs the scripts that need only files in this repository and compares what they produce with the bundled
# results. Runs in a temporary copy of the tracked files, so the working tree is never changed.
#
# Usage: bash ci/reproduce.sh [r|python|all]     (default all)
#   r       needs R with data.table and limma (CI: Ubuntu 24.04 packages, R 4.3.3, limma 3.58.1)
#   python  needs Python with numpy, pandas and matplotlib
#
# Printouts (.out) are compared after normalising whitespace and dropping data.table's column-type row
# (<char> <int> ..., printed by data.table >= 1.15 only) and the "R version" line; tables are compared byte for byte.
# To add a folder: add check_out / same_file / run_ok lines below and set ci=reproduce (or smoke) in catalog.csv.
set -uo pipefail
GROUP=${1:-all}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
(cd "$ROOT" && git ls-files -z | xargs -0 cp --parents -t "$WORK")
fail=0

norm() { sed -E '/^[[:space:]]*(<[a-z0-9]+>[[:space:]]*)+$/d; /^\[1\] "R version/d; s/[[:space:]]+/ /g; s/^ //; s/ $//' "$1"; }

# check_out <folder> <command> <bundled printout>: run the command in the folder, compare its printout
check_out() {
  local dir=$1 cmd=$2 ref=$3 log="$WORK/.ci.$RANDOM.out"
  if ! (cd "$WORK/$dir" && eval "$cmd") > "$log" 2>&1; then
    echo "FAIL $dir: '$cmd' exited with an error"; tail -20 "$log"; fail=1; return
  fi
  if diff <(norm "$ROOT/$dir/$ref") <(norm "$log") > "$log.diff"; then
    echo "ok   $dir/$ref"
  else
    echo "FAIL $dir/$ref: printout differs (< bundled, > re-run)"; head -40 "$log.diff"; fail=1
  fi
}

# same_file <path>: a file rewritten by a script must be byte for byte the bundled one
same_file() {
  if cmp -s "$ROOT/$1" "$WORK/$1"; then echo "ok   $1"; else echo "FAIL $1 differs from the bundled file"; fail=1; fi
}

# run_ok <folder> <command>: the command must finish without an error (outputs such as figures are not compared)
run_ok() {
  local dir=$1 cmd=$2 log="$WORK/.ci.$RANDOM.log"
  if (cd "$WORK/$dir" && eval "$cmd") > "$log" 2>&1; then
    echo "ok   $dir: $cmd"
  else
    echo "FAIL $dir: '$cmd' exited with an error"; tail -20 "$log"; fail=1
  fi
}

if [ "$GROUP" = r ] || [ "$GROUP" = all ]; then
  echo "== R syntax"
  if (cd "$ROOT" && git ls-files -z '*.R' | xargs -0 Rscript -e 'for (f in commandArgs(TRUE)) invisible(parse(f))'); then
    echo "ok   all .R files parse"
  else
    fail=1
  fi
  echo "== R re-runs"
  Rscript -e 'cat(R.version.string, "| data.table", format(packageVersion("data.table")), "| limma", format(packageVersion("limma")), "\n")'
  check_out diann-hye "Rscript score_diann.R" score_diann.out
  same_file diann-hye/diann_scores.csv
  check_out diann-2.7-hye "Rscript score270.R" score270.out
  same_file diann-2.7-hye/diann_scores_270.csv
  check_out diann-2.7-hye "Rscript common.R" common.out
  check_out diann-2.7-hye "Rscript filt.R" filt.out
  run_ok diann-2.7-hye "Rscript export_fc.R"
  same_file diann-2.7-hye/fc_both_versions.csv
  check_out de-pipeline-check/diann270 "Rscript diann_limma.R" diann_limma.out
  same_file de-pipeline-check/diann270/limma_results_filtered.csv
  check_out limma-3v3-sim "Rscript three_vs_three.R" output.txt
fi

if [ "$GROUP" = python ] || [ "$GROUP" = all ]; then
  echo "== Python re-runs"
  python3 -c 'import sys, numpy, pandas, matplotlib; print("Python", sys.version.split()[0], "| numpy", numpy.__version__, "| pandas", pandas.__version__, "| matplotlib", matplotlib.__version__)'
  run_ok string-hubs "python3 summarize_api.py"
  same_file string-hubs/api_summary_extra.json
  run_ok docs "MPLBACKEND=Agg python3 make_summary_figure.py"
  run_ok preprint-hye "MPLBACKEND=Agg python3 figures.py"
fi

if [ "$fail" -ne 0 ]; then echo; echo "Some checks failed."; exit 1; fi
echo; echo "All checks passed."
