#!/usr/bin/env bash
# Re-runs this folder's analysis from the raw inputs and rewrites its outputs next to this script.
# Usage: ./run.sh <input dir>      (the inputs are not in the repository; see README "Data")
set -euo pipefail
cd "$(dirname "$0")"
IN=${1:?usage: ./run.sh <input dir>}

# One line per step; each step's printout goes to <script>.out next to the script.
Rscript score.R "$IN" > score.out 2>&1
# python3 figure.py > figure.out 2>&1

# For long external runs whose exit status you log, capture it in the same command (set -e would otherwise stop the
# script before the next line, and $(date) resets $?, so "exit $?" written after it always prints 0):
#   rc=0; tool --args > run.log 2>&1 || rc=$?; echo "end $(date) exit $rc"
