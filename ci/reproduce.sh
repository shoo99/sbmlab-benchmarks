#!/usr/bin/env bash
# Re-runs the scripts that need only files in this repository and compares what they produce with the bundled
# results. Runs in a temporary copy of the tracked files, so the working tree is never changed.
#
# Usage: bash ci/reproduce.sh [r|python|all|pipelines|microbiome|doublets]     (default all = r + python)
#   r          needs R with data.table, limma, sva, bladderbatch, survival, maxstat, edgeR and DESeq2 (CI: Ubuntu 24.04
#              packages, R 4.3.3, limma 3.58.1)
#   python     needs Python with numpy, pandas and matplotlib
#   pipelines  end-to-end tests on synthetic data; needs salmon, kallisto, STAR, subread, RSEM, R with data.table
#   microbiome re-runs microbiome-da-methods in its pinned conda environment; set DA_RSCRIPT to that env's Rscript
#   doublets   the scrna-doublets pipeline on synthetic data in its pinned conda environment; set DBL_RSCRIPT and
#              DBL_PYTHON to that env's Rscript and python3
#
# Printouts (.out) are compared after normalising whitespace and dropping data.table's column-type row
# (<char> <int> ..., printed by data.table >= 1.15 only) and the "R version" line; tables are compared byte for byte.
# To add a folder: add check_out / same_file / run_ok lines below and set ci=reproduce (or smoke) in catalog.csv.
set -uo pipefail
GROUP=${1:-all}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
if ! (cd "$ROOT" && git ls-files -z) > "$WORK/.tracked" || [ ! -s "$WORK/.tracked" ]; then
  echo "Cannot list the tracked files with git ls-files (see the git message above)."
  echo "In a container, add the checkout as a safe directory: git config --global --add safe.directory \"\$PWD\""
  exit 1
fi
(cd "$ROOT" && xargs -0 cp --parents -t "$WORK" < "$WORK/.tracked") || { echo "Copying the tracked files failed."; exit 1; }
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
  if (cd "$WORK" && tr '\0' '\n' < .tracked | grep '\.R$' | tr '\n' '\0' | xargs -0 Rscript -e 'for (f in commandArgs(TRUE)) invisible(parse(f))'); then
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
  check_out batch-correction-unbalanced "Rscript batch_study.R" batch_study.out
  same_file batch-correction-unbalanced/scores.csv
  check_out random-signatures-tcga "OPENBLAS_NUM_THREADS=1 Rscript cutpoint_sim.R" cutpoint_sim.out
  same_file random-signatures-tcga/sim_scores.csv
  check_out random-signatures-tcga "Rscript score.R" score.out
  same_file random-signatures-tcga/scores.csv
  run_ok random-signatures-tcga "Rscript test/test_synthetic.R"     # tcga_signatures.R on synthetic data, known answers
  check_out scrna-pseudobulk-de "OPENBLAS_NUM_THREADS=1 Rscript sim_study.R" sim_study.out
  same_file scrna-pseudobulk-de/sim_scores.csv
  check_out scrna-pseudobulk-de "Rscript score.R" score.out
  same_file scrna-pseudobulk-de/scores.csv
  run_ok scrna-pseudobulk-de "Rscript test/test_kang_format.R"     # kang_study.R on synthetic data in the GEO layout
  run_ok common "Rscript tests/test_known_ratio.R"
fi

if [ "$GROUP" = python ] || [ "$GROUP" = all ]; then
  echo "== Python re-runs"
  python3 -c 'import sys, numpy, pandas, matplotlib; print("Python", sys.version.split()[0], "| numpy", numpy.__version__, "| pandas", pandas.__version__, "| matplotlib", matplotlib.__version__)'
  run_ok string-hubs "python3 summarize_api.py"
  same_file string-hubs/api_summary_extra.json
  run_ok docs "MPLBACKEND=Agg python3 make_summary_figure.py"
  run_ok preprint-hye "MPLBACKEND=Agg python3 figures.py"
  run_ok batch-correction-unbalanced "MPLBACKEND=Agg python3 figure.py"
  run_ok microbiome-da-methods "MPLBACKEND=Agg python3 figure.py"
  run_ok random-signatures-tcga "MPLBACKEND=Agg python3 figure.py"
  run_ok scrna-pseudobulk-de "MPLBACKEND=Agg python3 figure.py"
  check_out db-api-monitor "python3 summarize.py" summarize.out
  same_file db-api-monitor/scores.csv
  run_ok db-api-monitor "python3 test/test_monitor.py"       # monitor.py against a local server
  check_out metagenome-profilers-mock "python3 score.py" score.out
  same_file metagenome-profilers-mock/scores.csv
  run_ok metagenome-profilers-mock "MPLBACKEND=Agg python3 figure.py"
  check_out variant-calling-giab "python3 score.py" score.out
  same_file variant-calling-giab/scores.csv
  run_ok variant-calling-giab "MPLBACKEND=Agg python3 figure.py"
  run_ok common "python3 tests/test_known_ratio.py"
  run_ok entrapment-fdr-hye "python3 test/test_synthetic.py"    # entrapment database, extraction and FDP estimates
  run_ok docs "python3 make_summaries.py"            # index pages must match catalog.csv and every scores.csv
  same_file docs/catalog.md
  same_file docs/scores_summary.csv
  same_file docs/scores_summary.md
fi

if [ "$GROUP" = microbiome ]; then
  echo "== Microbiome DA study in its pinned conda environment (microbiome-da-methods/conda-linux-64.lock)"
  RS=${DA_RSCRIPT:-Rscript}
  "$RS" -e 'cat(R.version.string, "| ANCOMBC", format(packageVersion("ANCOMBC")), "| ALDEx2", format(packageVersion("ALDEx2")), "\n")'
  check_out microbiome-da-methods "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \"$RS\" da_study.R" da_study.out
  same_file microbiome-da-methods/scores.csv
fi

if [ "$GROUP" = doublets ]; then
  echo "== Doublet pipeline on synthetic data in its pinned conda environment (scrna-doublets/conda-linux-64.lock)"
  run_ok scrna-doublets "bash test/test_pipeline.sh"
fi

if [ "$GROUP" = pipelines ]; then
  echo "== Pipelines on synthetic data with a known answer"
  run_ok rnaseq-quant-seqc "bash test/test_pipeline.sh 4"
fi

if [ "$fail" -ne 0 ]; then echo; echo "Some checks failed."; exit 1; fi
echo; echo "All checks passed."
