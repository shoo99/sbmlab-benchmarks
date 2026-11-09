"""Figure for the Kang et al. 2018 run, drawn from kang_metrics.csv (written by kang_study.R).
Usage: python3 figure_kang.py   -> scrna_pseudobulk_kang.png next to this script (needs matplotlib)
  a  null (4 vs 4 donor splits, no real difference): share of genes with p < 0.05, per cell type
  b  planted 2-fold changes in 10% of genes: false discovery proportion of BH 5% calls
  c  planted: power of BH 5% calls
Only the default cell-level test (Wilcoxon) and one pseudobulk test (edgeR QL) are drawn; the t-test, DESeq2 and
limma-voom rows are in score.out and scores.csv.
"""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
SERIES = {"Wilcoxon (cells)": "#eb6834", "edgeR QL (pseudobulk)": "#2a78d6"}
LEGEND = {"Wilcoxon (cells)": "Wilcoxon on cells (Seurat default)", "edgeR QL (pseudobulk)": "edgeR QL on pseudobulk"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
TYPES = ["CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", "FCGR3A+ Monocytes"]
LABEL = {"CD4 T cells": "CD4 T", "CD14+ Monocytes": "CD14+\nmono", "B cells": "B", "CD8 T cells": "CD8 T", "FCGR3A+ Monocytes": "FCGR3A+\nmono"}

acc = defaultdict(list)
cells = {}
for r in csv.DictReader(open(HERE / "kang_metrics.csv")):
    cells[r["cell_type"]] = int(r["cells"])
    for k in ("p05", "fdp", "power"):
        if r[k] not in ("NA", ""):
            acc[(r["cell_type"], r["scenario"], r["method"], k)].append(float(r[k]))
mean = lambda ct, sc, m, k: sum(acc[(ct, sc, m, k)]) / len(acc[(ct, sc, m, k)])


def panel(ax, scenario, metric, title, ylabel, ref=None):
    w = 0.38
    for j, (m, colour) in enumerate(SERIES.items()):
        xs = [i + (j - 0.5) * (w + 0.03) for i in range(len(TYPES))]
        ys = [mean(ct, scenario, m, metric) for ct in TYPES]
        ax.bar(xs, ys, w, color=colour, label=LEGEND[m])
        for x, y in zip(xs, ys):
            ax.text(x, y + 0.015, f"{y:.2f}", ha="center", va="bottom", fontsize=6.5, color=MUTED)
    if ref is not None:
        ax.axhline(ref, color=INK, lw=0.8, ls="--")
    ax.set_xticks(range(len(TYPES)), [f"{LABEL[ct]}\n({cells[ct]:,} cells)" for ct in TYPES], fontsize=7)
    ax.set_title(title, fontsize=9, loc="left", color=INK)
    ax.set_ylabel(ylabel, fontsize=8, color=MUTED)
    ax.tick_params(labelsize=7, colors=MUTED)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)


fig, axs = plt.subplots(1, 3, figsize=(13, 4.1))
panel(axs[0], "null", "p05", "a  Same condition, donors split 4 vs 4 (35 splits)",
      "share of genes with p < 0.05  (dashed: 0.05 expected)", 0.05)
panel(axs[1], "planted", "fdp", "b  10% of genes changed 2-fold: false share of calls",
      "FDP of BH 5% calls  (dashed: 0.05 target)", 0.05)
panel(axs[2], "planted", "power", "c  10% of genes changed 2-fold: power", "changed genes found (BH 5%)")
axs[0].legend(fontsize=7, frameon=False, loc="upper right")
axs[0].set_ylim(0, 0.3)
axs[1].set_ylim(0, 0.7)
axs[2].set_ylim(0, 1.1)
fig.suptitle("Kang et al. 2018, unstimulated PBMC from 8 lupus patients (GSE96583 batch 2)", fontsize=9, x=0.01, ha="left", color=MUTED)
fig.tight_layout()
fig.savefig(HERE / "scrna_pseudobulk_kang.png", dpi=150)
print("wrote scrna_pseudobulk_kang.png")
