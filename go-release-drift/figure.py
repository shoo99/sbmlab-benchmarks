"""Figure for the GO release drift study, drawn from scores.csv and results/cholesterol_genes.tsv.
Usage: python3 figure.py   -> go_release_drift.png next to this script (needs matplotlib)
  a  significant terms per release for the same lists (solid: all evidence codes; dashed: without IEA)
  b  Jaccard overlap of the significant terms with those of the 2025-10-10 release
  c  ontology change alone (2016 annotations, 2016 vs 2025 ontology) and annotation change alone (2025 ontology,
     2016 vs 2025 annotations), against both together
  d  which detected cholesterol-biosynthesis proteins carry "cholesterol biosynthetic process" in each release
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = Path(__file__).resolve().parent
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
COLORS = {"planted lists (5)": "#1baf7a", "airway up": "#eb6834", "airway down": "#2a78d6"}
RELEASES = ["2016-01-01", "2020-01-01", "2023-01-01", "2025-10-10"]
EVID = {"all evidence codes": "-", "without IEA": "--"}
rows = list(csv.DictReader(open(HERE / "scores.csv")))
chol = list(csv.DictReader(open(HERE / "results" / "cholesterol_genes.tsv"), delimiter="\t"))


def val(version, setting, subset, metric):
    for r in rows:
        if r["version"] == version and r["setting"] == setting and r["subset"] == subset and r["metric"] == metric:
            return float(r["value"])
    return None


def style(ax, title):
    ax.set_title(title, fontsize=9, loc="left", color=INK)
    ax.tick_params(labelsize=7, colors=MUTED)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)


fig, axs = plt.subplots(2, 2, figsize=(12, 8))
ax = axs[0, 0]
for g, c in COLORS.items():
    for ev, ls in EVID.items():
        ys = [val(f"GO {r}", ev, g, "terms_found") for r in RELEASES]
        ax.plot(range(4), ys, ls=ls, color=c, marker="o", ms=4, lw=1.4, label=f"{g}, {ev}")
ax.set_xticks(range(4), RELEASES)
ax.set_ylabel("significant GO terms (FDR < 0.05)", fontsize=8, color=MUTED)
fr = [int(val(f"GO {r}", "all evidence codes", "random lists (20)", "false_terms")) for r in RELEASES]
ax.text(0.02, 0.62, "random lists (20 x 300 proteins): false terms " + ", ".join(f"{r[:4]} {n}" for r, n in zip(RELEASES, fr)),
        transform=ax.transAxes, fontsize=7, color=MUTED, va="top")
style(ax, "a  Terms found for unchanged lists")

ax = axs[0, 1]
for g, c in COLORS.items():
    for ev, ls in EVID.items():
        ys = [val(f"GO {r}", ev, g, "jaccard_vs_2025") for r in RELEASES[:-1]] + [1.0]
        ax.plot(range(4), ys, ls=ls, color=c, marker="o", ms=4, lw=1.4)
ax.set_xticks(range(4), RELEASES)
ax.set_ylim(0, 1.05)
ax.set_ylabel("Jaccard with the 2025-10-10 terms", fontsize=8, color=MUTED)
ax.legend(*axs[0, 0].get_legend_handles_labels(), fontsize=6.5, frameon=False, loc="lower right")
style(ax, "b  Overlap with the newest release")

ax = axs[1, 0]
parts = [("ontology alone", "#86b6ef"), ("annotations alone", "#2a78d6"), ("both (2016 vs 2025)", INK)]
hyb = "GO 2016-01-01 annotations, 2025-10-10 ontology"
for i, g in enumerate(COLORS):
    vals = [val(hyb, "all evidence codes", g, "jaccard_vs_2016"), val(hyb, "all evidence codes", g, "jaccard_vs_2025"),
            val("GO 2016-01-01", "all evidence codes", g, "jaccard_vs_2025")]
    for k, (v, (lab, c)) in enumerate(zip(vals, parts)):
        ax.bar(i + (k - 1) * 0.26, v, 0.24, color=c, label=lab if i == 0 else None)
        ax.text(i + (k - 1) * 0.26, v + 0.02, f"{v:.2f}", ha="center", fontsize=6.5, color=MUTED)
ax.set_xticks(range(3), list(COLORS))
ax.set_ylim(0, 1.05)
ax.set_ylabel("Jaccard of the significant terms", fontsize=8, color=MUTED)
ax.legend(fontsize=7, frameon=False, loc="upper right")
style(ax, "c  What changed: the ontology or the annotations (all evidence codes)")

ax = axs[1, 1]
genes = sorted({r["gene"] for r in chol})
cell = {(r["gene"], r["release"]): (r["GO0006695"], r["GO0016126"]) for r in chol if r["evidence"] == "all"}
for j, rel in enumerate(RELEASES):
    for i, g in enumerate(genes):
        c6695, c16126 = cell[(g, rel)]
        color = "#1baf7a" if c6695 == "1" else ("#a8dcc6" if c16126 == "1" else "#f1f0ec")
        ax.add_patch(plt.Rectangle((j, len(genes) - 1 - i), 0.92, 0.88, color=color))
ax.set_xlim(-0.05, 4)
ax.set_ylim(-0.1, len(genes))
ax.set_xticks([j + 0.46 for j in range(4)], RELEASES)
ax.set_yticks([len(genes) - 1 - i + 0.44 for i in range(len(genes))], genes, fontsize=6.5)
ax.legend(handles=[Patch(color="#1baf7a", label="cholesterol biosynthetic process"),
                   Patch(color="#a8dcc6", label="sterol biosynthetic process only"),
                   Patch(color="#f1f0ec", label="neither")], fontsize=6.5, frameon=False,
          loc="upper left", bbox_to_anchor=(1.0, 1.0))
style(ax, "d  Planted proteins (Reactome) carrying the GO term, all evidence codes")
ax.grid(False)
fig.tight_layout()
fig.savefig(HERE / "go_release_drift.png", dpi=150)
print("wrote go_release_drift.png")
