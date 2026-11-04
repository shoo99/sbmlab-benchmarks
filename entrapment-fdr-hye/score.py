"""Entrapment estimates of the false discovery proportion from counts.csv (written by extract.py); issue #5.

With N_T target and N_E entrapment identifications at a threshold, and r = entrapment / target database size
(peptides for precursor and peptide level, proteins for protein level), following Wen et al. (2025, Nat Methods):
  fdp_lower      N_E / (N_T + N_E)                 a lower bound: every entrapment hit is false, some targets are too
  fdp_combined   N_E (1 + 1/r) / (N_T + N_E)       the combined estimate: false targets are expected at N_E / r
A tool controls its FDR if fdp_combined stays at or below the nominal level; fdp_lower above it shows a failure.
Also: n_target, n_entrapment, and share_targets_lost = 1 - N_T / (targets of the target-only search) when given.
Usage: python3 score.py [counts.csv] [out_dir]   (defaults: this folder; standard library only)
Writes scores.csv into out_dir and prints a table.
"""
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
counts = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "counts.csv"
out = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE
FOLDER = "entrapment-fdr-hye"
COLUMNS = ["folder", "dataset", "tool", "version", "setting", "subset", "metric", "value", "n"]
DATASET = "PXD028735 HYE (6 runs), target + entrapment database"

rows = []
print(f"{'search':<45} {'level':<9} {'q':<8} {'targets':>8} {'entrap':>7} {'r':>6} {'fdp_lower':>9} {'fdp_comb':>9} {'lost':>6}")
for c in csv.DictReader(open(counts)):
    nt, ne, r = int(c["n_target"]), int(c["n_entrapment"]), float(c["r"])
    n = nt + ne
    m = {"n_target": nt, "n_entrapment": ne,
         "fdp_lower": round(ne / n, 5) if n else 0.0,
         "fdp_combined": round(ne * (1 + 1 / r) / n, 5) if n else 0.0}
    if c["n_baseline_target"]:
        m["share_targets_lost"] = round(1 - nt / int(c["n_baseline_target"]), 4)
    setting = f"MBR {c['mbr']}; q <= {c['threshold']}" if c["threshold"] != "reported" else f"MBR {c['mbr']}; as reported"
    for metric, v in m.items():
        rows.append([FOLDER, f"{DATASET}, r = {r:g}", c["search"], c["version"] or c["engine"], setting, c["level"],
                     metric, v, n])
    lost = f"{m['share_targets_lost']:.3f}" if "share_targets_lost" in m else "-"
    print(f"{c['search']:<45} {c['level']:<9} {c['threshold']:<8} {nt:>8} {ne:>7} {r:>6.3g} {m['fdp_lower']:>9.4f} "
          f"{m['fdp_combined']:>9.4f} {lost:>6}")
out.mkdir(parents=True, exist_ok=True)
with open(out / "scores.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(COLUMNS)
    w.writerows(rows)
print(f"wrote scores.csv ({len(rows)} rows)")
