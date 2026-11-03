"""Check that score_seqc.R recovers the simulated answer for every quantifier (standard library only).
Usage: python3 check_synthetic.py <scores.csv> <truth.json>

Thresholds are deliberately loose: the synthetic data are small (16 spikes, ~120 genes, 40k read pairs per
sample, replicate noise), so this checks that the pipeline and the scoring work, not how good each tool is.
"""
import csv
import json
import sys
from collections import defaultdict

scores, truth = sys.argv[1], json.load(open(sys.argv[2]))
by_tool = defaultdict(dict)
ercc = defaultdict(lambda: [0.0, 0])
for r in csv.DictReader(open(scores)):
    if r["metric"] == "within_0.5_centred":
        ercc[r["tool"]][0] += float(r["value"]) * int(r["n"])
        ercc[r["tool"]][1] += int(r["n"])
    else:
        by_tool[r["tool"]][r["metric"]] = float(r["value"])

fails = []
tools = {"salmon", "kallisto", "RSEM", "featureCounts"}
if set(by_tool) != tools:
    fails.append(f"tools scored: {sorted(by_tool)}; expected {sorted(tools)}")
for tool in sorted(by_tool):
    m = by_tool[tool]
    e = ercc[tool][0] / max(ercc[tool][1], 1)
    checks = [
        (f"ERCC within +-0.5 after centring {e:.2f} >= 0.80", e >= 0.80),
        (f"titration order consistent {m['titration_order_consistent']:.2f} >= 0.75", m["titration_order_consistent"] >= 0.75),
        (f"fitted w {m['mixture_w']:.3f} within 15% of {truth['w']}", abs(m["mixture_w"] / truth["w"] - 1) <= 0.15),
        (f"median |residual| C {m['mixture_median_abs_resid_C']:.3f}, D {m['mixture_median_abs_resid_D']:.3f} <= 0.30",
         max(m["mixture_median_abs_resid_C"], m["mixture_median_abs_resid_D"]) <= 0.30),
    ]
    for desc, ok in checks:
        print(f"{'ok  ' if ok else 'FAIL'} {tool}: {desc}")
        if not ok:
            fails.append(f"{tool}: {desc}")
if fails:
    print(f"\n{len(fails)} check(s) failed")
    sys.exit(1)
print("\nsynthetic answer recovered by all four quantifiers")
