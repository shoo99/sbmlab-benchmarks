"""Figure for the Multiome pairing benchmark, drawn from scores.csv (written by score.py).
Usage: python3 figure.py   -> multiome_integration_pairing.png next to this script (needs matplotlib)
  a  FOSCTTM per method and cell type (lower is better; 0.5 is random)
  b  label-transfer accuracy per method and cell type
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
PALETTE = ["#86b6ef", "#2a78d6", "#1baf7a", "#eb6834"]
rows = list(csv.DictReader(open(HERE / "scores.csv")))
methods = []
for r in rows:
    if r["tool"] not in methods:
        methods.append(r["tool"])
types = []
for r in rows:
    if r["subset"] != "all nuclei" and r["subset"] not in types:
        types.append(r["subset"])
cols = ["all nuclei"] + types


def val(m, subset, metric):
    for r in rows:
        if r["tool"] == m and r["subset"] == subset and r["metric"] == metric:
            return float(r["value"])
    return None


fig, axs = plt.subplots(2, 1, figsize=(12, 7.4))
w = 0.8 / len(methods)
for ax, metric, title, ylabel in ((axs[0], "foscttm", "a  FOSCTTM: share of the other modality closer than the true match",
                                   "FOSCTTM (lower is better; 0.5 = random)"),
                                  (axs[1], "label_accuracy", "b  Cell type transferred from RNA to the ATAC profile",
                                   "share of nuclei given their own RNA type")):
    for k, m in enumerate(methods):
        ys = [val(m, c, metric) or 0 for c in cols]
        ax.bar([i + (k - (len(methods) - 1) / 2) * w for i in range(len(cols))], ys, w * 0.92,
               color=PALETTE[k % len(PALETTE)], label=m)
    n = {c: next((r["n"] for r in rows if r["subset"] == c and r["metric"] == metric), "") for c in cols}
    ax.set_xticks(range(len(cols)), [f"{c}\n({n[c]})" for c in cols], fontsize=7)
    ax.set_title(title, fontsize=9, loc="left", color=INK)
    ax.set_ylabel(ylabel, fontsize=8, color=MUTED)
    ax.tick_params(labelsize=7, colors=MUTED)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)
axs[0].axhline(0.5, color=INK, lw=0.8, ls="--")
axs[1].set_ylim(0, 1.05)
axs[0].legend(fontsize=7, frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(HERE / "multiome_integration_pairing.png", dpi=150)
print("wrote multiome_integration_pairing.png")
