"""Convert Thermo Fisher's ERCC concentration table (cms_095046.txt) into ercc_truth.tsv for score_seqc.R.
Columns are found by their header text ("ERCC ID", "subgroup", "log2(Mix 1/Mix 2)"), so column order does not matter.
Usage: python3 ercc_truth.py <concentration table> <ercc_truth.tsv>
"""
import csv
import math
import sys

src, dst = sys.argv[1], sys.argv[2]
with open(src, newline="", encoding="utf-8-sig", errors="replace") as fh:
    rows = [[c.strip() for c in r] for r in csv.reader(fh, delimiter="\t") if any(c.strip() for c in r)]
header = [h.lower() for h in rows[0]]


def col(*words):
    hits = [i for i, h in enumerate(header) if all(w in h for w in words)]
    if len(hits) != 1:
        sys.exit(f"{src}: expected one column with {words} in its header, found {len(hits)}: {rows[0]}")
    return hits[0]


i_id, i_sub, i_log2 = col("ercc", "id"), col("subgroup"), col("log2")
out = [(r[i_id], r[i_sub], float(r[i_log2])) for r in rows[1:] if r[i_id].startswith("ERCC-")]
expected = {"A": 2.0, "B": 0.0, "C": math.log2(1 / 1.5), "D": -1.0}
if len(out) != 92:
    sys.exit(f"{src}: {len(out)} ERCC rows, expected 92")
for e, sub, v in out:
    if sub not in expected or abs(v - expected[sub]) > 0.01:
        sys.exit(f"{src}: {e} subgroup {sub} log2 ratio {v} does not match the ERCC Mix 1 / Mix 2 design")
with open(dst, "w") as fh:
    fh.write("ercc_id\tsubgroup\tlog2_mix1_over_mix2\n")
    fh.writelines(f"{e}\t{sub}\t{v:.6f}\n" for e, sub, v in out)
print(f"{dst}: 92 ERCC spikes, subgroups " + ", ".join(f"{s} ({sum(1 for _, x, _ in out if x == s)})" for s in "ABCD"))
