"""Figure for the random-signature study, drawn from scores.csv (written by score.R).
Usage: python3 figure.py   -> random_signatures_tcga.png next to this script (needs matplotlib)
  a  simulation, marker unrelated to survival: share of repetitions with p < 0.05 by analysis and sample size
  b  TCGA, random 50-gene sets with permuted survival labels (the truth is no association)
  c  TCGA, random 50-gene sets with the observed labels: all genes vs proliferation-free genes
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
SIM = ["Cox, continuous score", "median split", "optimal cut-point, naive p", "optimal cut-point, corrected p"]
TCGA = SIM + ["Cox, adjusted for proliferation", "Cox, adjusted for proliferation and stage"]
SHORT = {"Cox, continuous score": "Cox\ncontinuous", "median split": "median\nsplit",
         "optimal cut-point, naive p": "optimal cut\nnaive p", "optimal cut-point, corrected p": "optimal cut\ncorrected p",
         "Cox, adjusted for proliferation": "Cox adj.\nprolif.", "Cox, adjusted for proliferation and stage": "Cox adj.\nprolif.+stage"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
N_STEPS = {100: "#86b6ef", 300: "#2a78d6", 1000: "#104281"}           # ordinal blue ramp (light -> dark)
FAMILY = {"random 50-gene sets": "#2a78d6", "random 50-gene sets, proliferation-free": "#eb6834"}
SHAPES = ["o", "s", "^", "D"]

rows = list(csv.DictReader(open(HERE / "scores.csv")))
cohorts = []
for r in rows:
    if r["dataset"].startswith("TCGA "):
        c = r["dataset"].split()[1].rstrip(",")
        if c not in cohorts:
            cohorts.append(c)


def value(dataset_start, tool, setting_start, subset, metric="share_p_below_0.05"):
    for r in rows:
        if (r["dataset"].startswith(dataset_start) and r["tool"] == tool and r["setting"].startswith(setting_start)
                and r["subset"] == subset and r["metric"] == metric):
            return float(r["value"])
    return None


def style(ax, labels, title, ylabel):
    ax.axhline(0.05, color=INK, lw=0.8, ls="--")
    ax.set_xticks(range(len(labels)), [SHORT[a] for a in labels], fontsize=7)
    ax.set_title(title, fontsize=9, loc="left", color=INK)
    ax.set_ylabel(ylabel + "  (dashed: 0.05 expected)", fontsize=8, color=MUTED)
    ax.tick_params(labelsize=7, colors=MUTED)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)


fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), gridspec_kw={"width_ratios": [4, 6, 6]})
ax = axs[0]
w = 0.26
for k, (n, colour) in enumerate(N_STEPS.items()):
    ys = [value(f"simulation: n = {n},", a, "null", "simulated marker") for a in SIM]
    ax.bar([i + (k - 1) * (w + 0.02) for i in range(len(SIM))], ys, w, color=colour, label=f"n = {n}")
style(ax, SIM, "a  Simulation, marker unrelated to survival", "share of repetitions with p < 0.05")
ax.legend(fontsize=7, frameon=False, loc="upper left")

off = [-0.075, -0.025, 0.025, 0.075]
ax = axs[1]
for j, c in enumerate(cohorts):
    ys = [value(f"TCGA {c},", a, "permuted", "random 50-gene sets") for a in TCGA]
    ax.scatter([i + 3 * off[j] for i in range(len(TCGA))], ys, s=36, marker=SHAPES[j], color=MUTED,
               edgecolor="white", linewidth=0.8, label=c, zorder=3)
style(ax, TCGA, "b  TCGA, random 50-gene sets, permuted survival labels", "share of sets with p < 0.05")
ax.legend(fontsize=7, frameon=False, loc="upper right", title="cohort (PFI)", title_fontsize=7)

ax = axs[2]
for fam, colour in FAMILY.items():
    shift = -0.2 if fam == "random 50-gene sets" else 0.2
    for j, c in enumerate(cohorts):
        ys = [value(f"TCGA {c},", a, "observed", fam) for a in TCGA]
        ax.scatter([i + shift + off[j] for i in range(len(TCGA))], ys, s=36, marker=SHAPES[j], color=colour,
                   edgecolor="white", linewidth=0.8, zorder=3, label=c if fam == "random 50-gene sets" else None)
style(ax, TCGA, "c  TCGA, random 50-gene sets, observed survival labels", "share of sets with p < 0.05")
handles = [plt.Line2D([], [], ls="", marker="o", ms=6, color=col, label=("all genes" if "free" not in f else "proliferation-free genes"))
           for f, col in FAMILY.items()]
handles += [plt.Line2D([], [], ls="", marker=SHAPES[j], ms=6, color=MUTED, label=c) for j, c in enumerate(cohorts)]
ax.legend(handles=handles, fontsize=7, frameon=False, loc="upper right")
ymax = max(a.get_ylim()[1] for a in axs)
for a in axs:
    a.set_ylim(0, ymax)
fig.tight_layout()
fig.savefig(HERE / "random_signatures_tcga.png", dpi=150)
print("wrote random_signatures_tcga.png")
