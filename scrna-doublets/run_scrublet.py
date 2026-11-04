"""Runs Scrublet on one prepared data set (from prepare.R) and writes scores_scrublet.tsv next to it:
barcode, tool, score (higher = more doublet-like), call (Scrublet's own threshold; NA if it finds none).

Scrublet(expected_doublet_rate = the expected rate in info.tsv, random_state = 0), scrub_doublets() with its
defaults (min_counts 3, min_cells 3, 85th percentile gene variability, 30 PCs).
Usage: python3 run_scrublet.py <work_dir>/<dataset>   (needs scrublet, scipy, numpy)
"""
import csv
import gzip
import sys
import time
from pathlib import Path

import numpy as np
import scipy
import scipy.io
import scipy.sparse
import scrublet

d = Path(sys.argv[1])
with gzip.open(d / "counts.mtx.gz", "rb") as fh:
    X = scipy.sparse.csc_matrix(scipy.io.mmread(fh).T)      # droplets x genes
barcodes = (d / "barcodes.tsv").read_text().split()
info = next(csv.DictReader(open(d / "info.tsv"), delimiter="\t"))
rate = float(info["expected_rate"])
t0 = time.time()
s = scrublet.Scrublet(X, expected_doublet_rate=rate, random_state=0)
scores, calls = s.scrub_doublets(verbose=False)
with open(d / "scores_scrublet.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["barcode", "tool", "score", "call"])
    for i, b in enumerate(barcodes):
        call = "NA" if calls is None else ("TRUE" if calls[i] else "FALSE")
        w.writerow([b, "Scrublet", f"{scores[i]:.5g}", call])
version = getattr(scrublet, "__version__", None)
if version is None:
    from importlib.metadata import version as pkg_version
    version = pkg_version("scrublet")
with open(d / "versions_scrublet.tsv", "w") as fh:
    fh.write(f"tool\tversion\tsetting\nScrublet\tscrublet {version}\texpected_doublet_rate {rate}\n")
print(f"{d.name}: {len(barcodes)} droplets, expected rate {rate:.3f}, threshold "
      f"{'none' if calls is None else f'{s.threshold_:.3f}'}; {time.time() - t0:.0f} s")
print(f"Python {sys.version.split()[0]} | scrublet {version} | numpy {np.__version__} | scipy {scipy.__version__}")
