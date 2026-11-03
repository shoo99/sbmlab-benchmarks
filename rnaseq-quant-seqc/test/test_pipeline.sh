#!/usr/bin/env bash
# End-to-end check on synthetic data with a known answer:
#   simulate/simulate.py -> run.sh -> gene_matrix.R -> score_seqc.R -> test/check_synthetic.py
# Needs salmon, kallisto, STAR, featureCounts (subread), RSEM, R with data.table, and python3.
# Runs in a temporary directory (about 2 minutes with 4 threads); nothing is written to the repository.
# Usage: bash test/test_pipeline.sh [threads]
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
T=${1:-4}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
python3 "$HERE/simulate/simulate.py" "$TMP/data" --pairs 40000 --seed 1
bash "$HERE/run.sh" "$TMP/data" "$T"
Rscript "$HERE/gene_matrix.R" "$TMP/data/work" "$TMP/out"
Rscript "$HERE/score_seqc.R" "$TMP/out" "$TMP/data/ref/ercc_truth.tsv" "$TMP/out" synthetic 3
python3 "$HERE/test/check_synthetic.py" "$TMP/out/scores.csv" "$TMP/data/truth.json"
