#!/bin/bash
# Full run: about 5 minutes on 16 threads. WGCNA goes into ./Rlib unless already installed.
set -e
W="$(cd "$(dirname "$0")" && pwd)"; cd "$W"; mkdir -p out logs
export R_LIBS_USER="${R_LIBS_USER:-$W/Rlib}"
Rscript -e 'quit(status = !requireNamespace("WGCNA", quietly = TRUE))' 2>/dev/null || Rscript scripts/install.R
Rscript scripts/01_prep.R "$W" > logs/01_prep.txt 2>&1
for v in A B C; do
  Rscript scripts/02_wgcna.R "$W" $v > logs/02_wgcna_$v.txt 2>&1
  Rscript scripts/03_cv.R "$W" $v > logs/03_cv_$v.txt 2>&1
done
echo done
