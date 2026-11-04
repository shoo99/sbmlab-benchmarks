"""Scores the unpaired integrations against the true Multiome pairing (results/<method>.tsv); issue #16.

Each results file has one row per nucleus: its RNA cell type (true_type), the type transferred to its ATAC profile
(predicted_type), and FOSCTTM in both directions (foscttm_atac: share of RNA cells closer to the ATAC profile than its
own RNA profile; foscttm_rna: the reverse). Metrics per method:
  foscttm                mean over both directions (0 = every true match is the nearest; 0.5 = random)
  foscttm_atac_to_rna, foscttm_rna_to_atac
  label_accuracy         share of ATAC profiles given their nucleus's RNA cell type
  errors_within_lineage  share of the label errors that stay within a lineage (T/NK, B, monocytes, dendritic cells),
                         e.g. CD4 T taken for CD8 T; the rest cross lineages
  per cell type: label accuracy and FOSCTTM
Usage: python3 score.py [results_dir] [out_dir]   (defaults: this folder; standard library only)
"""
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
res = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "results"
out = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE
FOLDER = "multiome-integration-pairing"
COLUMNS = ["folder", "dataset", "tool", "version", "setting", "subset", "metric", "value", "n"]
DATASET = "10x PBMC granulocyte-sorted 10k Multiome (RNA and ATAC of the same nuclei, pairing hidden)"
METHODS = {"naive_baseline": ("naive baseline (gene activity + joint PCA)", "R 4.3.3; irlba 2.3.5.1"),
           "seurat_cca": ("Seurat CCA (label transfer, co-embedding)", "Seurat 5.0.1"),
           "glue": ("GLUE", "scglue 0.5.2"), "liger": ("LIGER / UINMF", "rliger")}
LINEAGE = {"CD4 T": "T/NK", "CD8 T": "T/NK", "NK": "T/NK", "B": "B", "CD14 Mono": "monocytes",
           "CD16 Mono": "monocytes", "cDC": "dendritic cells", "pDC": "dendritic cells"}

rows = []
print(f"== {DATASET}\n")
print(f"   {'method':<44} {'FOSCTTM':>8} {'ATAC>RNA':>9} {'RNA>ATAC':>9} {'labels':>7} {'errors in lineage':>18}")
per_type = {}
for f in sorted(res.glob("*.tsv"), key=lambda p: list(METHODS).index(p.stem) if p.stem in METHODS else 99):
    label, version = METHODS.get(f.stem, (f.stem, ""))
    d = list(csv.DictReader(open(f), delimiter="\t"))
    n = len(d)
    fa = [float(r["foscttm_atac"]) for r in d]
    fr = [float(r["foscttm_rna"]) for r in d]
    ok = [r["true_type"] == r["predicted_type"] for r in d]
    err = [r for r in d if r["true_type"] != r["predicted_type"]]
    within = sum(LINEAGE.get(r["true_type"]) == LINEAGE.get(r["predicted_type"]) for r in err)
    m = {"foscttm": (sum(fa) + sum(fr)) / (2 * n), "foscttm_atac_to_rna": sum(fa) / n, "foscttm_rna_to_atac": sum(fr) / n,
         "label_accuracy": sum(ok) / n, "errors_within_lineage": within / len(err) if err else 0.0}
    for metric, v in m.items():
        rows.append([FOLDER, DATASET, label, version, "unpaired (pairing hidden)", "all nuclei", metric, round(v, 4), n])
    print(f"   {label:<44} {m['foscttm']:>8.4f} {m['foscttm_atac_to_rna']:>9.4f} {m['foscttm_rna_to_atac']:>9.4f} "
          f"{m['label_accuracy']:>7.3f} {m['errors_within_lineage']:>18.3f}")
    types = sorted({r["true_type"] for r in d})
    per_type[label] = {}
    for t in types:
        sub = [r for r in d if r["true_type"] == t]
        acc = sum(r["predicted_type"] == t for r in sub) / len(sub)
        fo = sum(float(r["foscttm_atac"]) + float(r["foscttm_rna"]) for r in sub) / (2 * len(sub))
        per_type[label][t] = (acc, fo, len(sub))
        rows.append([FOLDER, DATASET, label, version, "unpaired (pairing hidden)", t, "label_accuracy", round(acc, 4), len(sub)])
        rows.append([FOLDER, DATASET, label, version, "unpaired (pairing hidden)", t, "foscttm", round(fo, 4), len(sub)])
if per_type:
    labels = list(per_type)
    types = list(next(iter(per_type.values())))
    print(f"\n   per cell type: label accuracy / FOSCTTM")
    print(f"   {'type':<11} {'nuclei':>6}  " + "  ".join(f"{l[:26]:>26}" for l in labels))
    for t in types:
        cells = "  ".join(f"{per_type[l][t][0]:>12.3f} / {per_type[l][t][1]:<11.4f}" for l in labels)
        print(f"   {t:<11} {per_type[labels[0]][t][2]:>6}  {cells}")
out.mkdir(parents=True, exist_ok=True)
with open(out / "scores.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(COLUMNS)
    w.writerows(rows)
print(f"\nwrote scores.csv ({len(rows)} rows)")
