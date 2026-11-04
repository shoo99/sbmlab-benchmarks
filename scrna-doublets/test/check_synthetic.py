"""Checks of the synthetic doublet run (called by test/test_pipeline.sh; standard library only).
Usage: python3 check_synthetic.py <generator_truth.tsv> <bundle_dir> <scores.csv>
  truth bookkeeping  prepare.R's truth equals the generator's labels: demuxlet-style labels for the Kang-like lanes
                     (ambs -> unknown); every cross-species doublet called "doublet" and every pure droplet "singlet"
                     for the hgmm-like set
  detection          every tool's AUPRC is at least twice the doublet share (chance); scDblFinder and Scrublet find
                     at least half the visible doublets at the expected rate
"""
import csv
import gzip
import sys
from pathlib import Path

gen_path, bundle, scores_path = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
gen = {(r["dataset"], r["barcode"]): r for r in csv.DictReader(open(gen_path), delimiter="\t")}
with gzip.open(bundle / "cell_scores.tsv.gz", "rt") as fh:
    cells = list(csv.DictReader(fh, delimiter="\t"))
scores = list(csv.DictReader(open(scores_path)))
checks = []

for ds in ("kang_ctrl", "kang_stim"):
    rows = [c for c in cells if c["dataset"] == ds]
    want = {"doublet": "doublet", "singlet": "singlet", "ambs": "unknown"}
    bad = sum(1 for c in rows if want[gen[(ds, c["barcode"])]["label"]] != c["truth"])
    checks.append((f"{ds}: truth matches the demuxlet-style labels ({len(rows)} droplets)", bad == 0 and len(rows) == 1200))
rows = [c for c in cells if c["dataset"] == "hgmm_12k"]
cross = [c for c in rows if gen[("hgmm_12k", c["barcode"])]["visible"] == "1"]
pure = [c for c in rows if gen[("hgmm_12k", c["barcode"])]["doublet"] == "0"]
checks.append((f"hgmm_12k: all {len(cross)} cross-species doublets labelled doublet", all(c["truth"] == "doublet" for c in cross)))
checks.append((f"hgmm_12k: all {len(pure)} pure droplets labelled singlet", all(c["truth"] == "singlet" for c in pure)))

tools = sorted({r["tool"] for r in scores})
checks.append(("seven tools scored (scDblFinder x2, scds x3, DoubletFinder, Scrublet)", len(tools) == 7))
for r in scores:
    if r["metric"] != "auprc":
        continue
    share = next(float(s["value"]) for s in scores
                 if s["dataset"] == r["dataset"] and s["tool"] == r["tool"] and s["metric"] == "doublet_share")
    checks.append((f"{r['dataset'].split(':')[0]}: {r['tool']} AUPRC {float(r['value']):.2f} >= 2 x doublet share {share:.2f}",
                   float(r["value"]) >= 2 * share))
for r in scores:
    if r["metric"] == "recall_at_expected" and r["tool"] in ("scDblFinder (expected rate)", "Scrublet"):
        checks.append((f"{r['dataset'].split(':')[0]}: {r['tool']} recall at expected rate {float(r['value']):.2f} >= 0.5",
                       float(r["value"]) >= 0.5))

for name, ok in checks:
    print(("ok   " if ok else "FAIL ") + name)
if not all(ok for _, ok in checks):
    sys.exit(1)
print("all synthetic checks passed")
