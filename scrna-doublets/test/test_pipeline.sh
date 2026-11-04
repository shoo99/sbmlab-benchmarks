#!/usr/bin/env bash
# End-to-end known-answer test of the doublet pipeline on synthetic data in the download layouts (run by CI).
#   make_synthetic.R -> prepare.R -> run_r_tools.R + run_scrublet.py per data set -> collect.py -> score.py
#   -> check_synthetic.py (truth bookkeeping and detection well above chance)
# Usage: bash test/test_pipeline.sh
#   env: DBL_RSCRIPT, DBL_PYTHON (default Rscript, python3: an R with scDblFinder, scds, Seurat, DoubletFinder and a
#        Python with scrublet; CI builds both from conda-linux-64.lock); DF_SWEEP (default 1)
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1
RS=${DBL_RSCRIPT:-Rscript}; PY=${DBL_PYTHON:-python3}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$RS" test/make_synthetic.R "$T/data"
EXPECTED_RATE=0.1 "$RS" prepare.R "$T/data" "$T/work"     # the synthetic sets hold 8-10% doublets
for d in "$T"/work/*/; do
  "$RS" run_r_tools.R "$d"
  "$PY" run_scrublet.py "$d"
done
"$PY" collect.py "$T/work" "$T/bundle"
"$PY" score.py "$T/bundle" "$T/out"
"$PY" test/check_synthetic.py "$T/data/generator_truth.tsv" "$T/bundle" "$T/out/scores.csv"
