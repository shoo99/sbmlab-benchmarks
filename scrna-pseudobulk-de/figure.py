"""Figure for the pseudobulk study, drawn from scores.csv (written by score.R from sim_study.R's sim_scores.csv).
Usage: python3 figure.py   -> scrna_pseudobulk_de.png next to this script (needs matplotlib)
  a  null: share of genes with p < 0.05, with and without donor-to-donor variation
  b  planted 2-fold changes: false discovery proportion of BH 5% calls
  c  planted 2-fold changes: power of BH 5% calls
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
METHODS = ["Wilcoxon (cells)", "t-test (cells)", "edgeR QL (pseudobulk)", "DESeq2 (pseudobulk)", "limma-voom (pseudobulk)"]
LABEL = {m: m.replace(" (cells)", "\ncells").replace(" (pseudobulk)", "\npseudobulk") for m in METHODS}
SERIES = {"donor SD 0.25": "#2a78d6", "no donor effect": "#eb6834"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
rows = [r for r in csv.DictReader(open(HERE / "scores.csv")) if r["dataset"].startswith("simulation")]


def value(method, setting_start, donor, metric):
    for r in rows:
        if r["tool"] == method and r["setting"].startswith(setting_start) and r["setting"].endswith(donor) and r["metric"] == metric:
            return float(r["value"])
    return None


def panel(ax, setting, metric, title, ylabel, ref=None):
    w = 0.38
    for k, (donor, colour) in enumerate(SERIES.items()):
        ys = [value(m, setting, donor, metric) for m in METHODS]
        xs = [i + (k - 0.5) * (w + 0.03) for i in range(len(METHODS))]
        ax.bar(xs, ys, w, color=colour, label=donor)
    if ref is not None:
        ax.axhline(ref, color=INK, lw=0.8, ls="--")
    ax.set_xticks(range(len(METHODS)), [LABEL[m] for m in METHODS], fontsize=7)
    ax.set_title(title, fontsize=9, loc="left", color=INK)
    ax.set_ylabel(ylabel, fontsize=8, color=MUTED)
    ax.tick_params(labelsize=7, colors=MUTED)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)


fig, axs = plt.subplots(1, 3, figsize=(13, 3.9))
panel(axs[0], "null", "share_p_below_0.05", "a  No difference between donor groups",
      "share of genes with p < 0.05  (dashed: 0.05 expected)", 0.05)
panel(axs[1], "planted", "fdp_bh_0.05", "b  10% of genes 2-fold: false share of calls",
      "FDP of BH 5% calls  (dashed: 0.05 target)", 0.05)
panel(axs[2], "planted", "power_bh_0.05", "c  10% of genes 2-fold: power", "changed genes found (BH 5%)")
axs[0].legend(fontsize=7, frameon=False, loc="upper right", title="simulated donors", title_fontsize=7)
for a in axs[:2]:
    a.set_ylim(0, 1)
axs[2].set_ylim(0, 1.05)
fig.tight_layout()
fig.savefig(HERE / "scrna_pseudobulk_de.png", dpi=150)
print("wrote scrna_pseudobulk_de.png")
