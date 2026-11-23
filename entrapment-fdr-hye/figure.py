"""Figures for the entrapment benchmark, drawn from the bundled tables (needs matplotlib).
Usage: python3 figure.py   -> entrapment_fdr_hye.png and entrapment_fdr_hye_all.png next to this script
entrapment_fdr_hye.png (from scores.csv, Arabidopsis entrapment, the first six searches): one panel per level
(precursors, peptides, protein groups, per run, MBR transfers only); one row per search. The bar runs from the lower
bound (fdp_lower) to the combined estimate (fdp_combined, the dot) at the nominal 1% (DIA-NN q <= 0.01; MaxQuant and
FragPipe as reported, MaxQuant protein groups at Q-value <= 0.01). Filled dots: MBR on; open dots: MBR off. The dashed
line is the nominal 1%.
entrapment_fdr_hye_all.png (from scores.csv, scores_orbitrap_dia.csv, paired.csv, paired_orbitrap_dia.csv): one panel
per level (precursors, peptides, protein groups); one row per search (MaxQuant, FragPipe, DIA-NN SWATH, DIA-NN
Orbitrap DIA; MBR off and on), with two lines per row: Arabidopsis entrapment (upper, circle) and shuffled entrapment
(lower, square), each from the lower bound to the combined estimate, and the paired estimate of the shuffled design
(diamond; peptides and protein groups only), all at the nominal 1%.
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
SEARCHES = ["DIA-NN library-free, MBR off", "DIA-NN library-free, MBR on", "MaxQuant, MBR off", "MaxQuant, MBR on",
            "FragPipe LFQ-MBR, MBR off", "FragPipe LFQ-MBR, MBR on"]
LABEL = {"DIA-NN library-free, MBR off": "DIA-NN 2.3.2, MBR off", "DIA-NN library-free, MBR on": "DIA-NN 2.3.2, MBR on",
         "MaxQuant, MBR off": "MaxQuant 2.8.1, MBR off", "MaxQuant, MBR on": "MaxQuant 2.8.1, MBR on",
         "FragPipe LFQ-MBR, MBR off": "FragPipe 24.0, MBR off", "FragPipe LFQ-MBR, MBR on": "FragPipe 24.0, MBR on"}
COLOUR = {"DIA-NN": "#2a78d6", "MaxQuant": "#eb6834", "FragPipe": "#1baf7a"}
LEVELS = [("precursor", "Precursors"), ("peptide", "Peptides"), ("protein", "Protein groups"),
          ("precursor_run", "Per run (run x precursor)"), ("precursor_run_mbr", "MBR transfers only")]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"

rows = list(csv.DictReader(open(HERE / "scores.csv")))


def val(search, level, metric):
    for r in rows:
        nominal = r["setting"].endswith("q <= 0.01") or r["setting"].endswith("as reported")
        if r["tool"] == search and r["subset"] == level and r["metric"] == metric and nominal:
            return float(r["value"]) * 100
    return None


fig, axs = plt.subplots(1, len(LEVELS), figsize=(14, 3.6), sharey=True)
xmax = max(v for s in SEARCHES for lev, _ in LEVELS for v in [val(s, lev, "fdp_combined")] if v is not None) * 1.15
for ax, (level, title) in zip(axs, LEVELS):
    for i, s in enumerate(SEARCHES):
        y = len(SEARCHES) - 1 - i
        colour = COLOUR[s.split(" ")[0].split(",")[0]]
        lo, comb = val(s, level, "fdp_lower"), val(s, level, "fdp_combined")
        if comb is None:
            ax.text(0.05, y, "not marked" if s.startswith("DIA-NN") else "none (MBR off)", va="center", fontsize=7, color=MUTED)
            continue
        ax.plot([lo, comb], [y, y], color=colour, lw=2, solid_capstyle="round")
        on = s.endswith("MBR on")
        ax.plot([comb], [y], "o", ms=8, color=colour if on else "white", mec=colour, mew=2)
    ax.axvline(1.0, color=MUTED, lw=1, ls=(0, (4, 3)))
    ax.set_xlim(0, xmax)
    ax.set_title(title, fontsize=9, loc="left", color=INK)
    ax.set_xlabel("estimated FDP, %", fontsize=8, color=MUTED)
    ax.tick_params(labelsize=7, colors=MUTED)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)
axs[0].set_yticks(range(len(SEARCHES)), [LABEL[s] for s in reversed(SEARCHES)], fontsize=8)
axs[0].set_ylim(-0.6, len(SEARCHES) - 0.4)
fig.suptitle("Entrapment estimate of the false discovery proportion at a nominal 1% FDR (HYE, PXD028735). "
             "Bar: lower bound to combined estimate (dot); filled = MBR on, open = MBR off; dashed = 1%",
             fontsize=8.5, color=MUTED, x=0.01, ha="left")
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(HERE / "entrapment_fdr_hye.png", dpi=150)
print("wrote entrapment_fdr_hye.png")


# ---- entrapment_fdr_hye_all.png: both designs, the paired estimate and the DIA-NN Orbitrap DIA searches
def read(name):
    return list(csv.DictReader(open(HERE / name)))


all_scores = read("scores.csv") + read("scores_orbitrap_dia.csv")
all_paired = read("paired.csv") + read("paired_orbitrap_dia.csv")


def at1(search, level, metric):
    for r in all_scores:
        nominal = r["setting"].endswith("q <= 0.01") or r["setting"].endswith("as reported")
        if r["tool"] == search and r["subset"] == level and r["metric"] == metric and nominal:
            return float(r["value"]) * 100
    return None


def paired(search, level):
    for r in all_paired:
        if r["search"] == search and r["level"] == level:
            return float(r["fdp_paired"]) * 100
    return None


ALL = [("MaxQuant", "MaxQuant 2.8.1"), ("FragPipe LFQ-MBR", "FragPipe 24.0"),
       ("DIA-NN library-free", "DIA-NN 2.3.2 SWATH"), ("DIA-NN Orbitrap DIA", "DIA-NN 2.3.2 Orbitrap DIA")]
ALL_COLOUR = {"MaxQuant": "#eb6834", "FragPipe LFQ-MBR": "#1baf7a", "DIA-NN library-free": "#2a78d6",
              "DIA-NN Orbitrap DIA": "#eda100"}
rows_all = [(f"{name}, MBR {m}", f"{label}, MBR {m}", ALL_COLOUR[name]) for name, label in ALL for m in ("off", "on")]
SHUF = " (shuffled entrapment)"
LEVELS_ALL = [("precursor", "Precursors"), ("peptide", "Peptides"), ("protein", "Protein groups")]
values = [v for s, _, _ in rows_all for lev, _ in LEVELS_ALL
          for v in (at1(s, lev, "fdp_combined"), at1(s + SHUF, lev, "fdp_combined"), paired(s + SHUF, lev)) if v]
xmax_all = max(values) * 1.12
fig, axs = plt.subplots(1, len(LEVELS_ALL), figsize=(12, 6.2), sharey=True)
for ax, (level, title) in zip(axs, LEVELS_ALL):
    for i, (s, _, colour) in enumerate(rows_all):
        y = len(rows_all) - 1 - i
        for dy, name, marker in ((0.18, s, "o"), (-0.18, s + SHUF, "s")):
            lo, comb = at1(name, level, "fdp_lower"), at1(name, level, "fdp_combined")
            if comb is None:
                continue
            ax.plot([lo, comb], [y + dy, y + dy], color=colour, lw=2, solid_capstyle="round")
            ax.plot([comb], [y + dy], marker, ms=7, color=colour, mec="white", mew=1.5)
        pr = paired(s + SHUF, level)
        if pr is not None:
            ax.plot([pr], [y - 0.18], "D", ms=6, color=INK, mec="white", mew=1.2)
    for y in range(len(rows_all) - 1):
        if y % 2 == 1:
            ax.axhline(y + 0.5, color=GRID, lw=0.8)
    ax.axvline(1.0, color=MUTED, lw=1, ls=(0, (4, 3)))
    ax.set_xlim(0, xmax_all)
    ax.set_title(title, fontsize=9, loc="left", color=INK)
    ax.set_xlabel("estimated FDP, %", fontsize=8, color=MUTED)
    ax.tick_params(labelsize=7, colors=MUTED)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)
axs[0].set_yticks(range(len(rows_all)), [label for _, label, _ in reversed(rows_all)], fontsize=8)
axs[0].set_ylim(-0.6, len(rows_all) - 0.4)
handles = [plt.Line2D([], [], color=MUTED, lw=2, marker="o", ms=7, mec="white",
                      label="Arabidopsis entrapment: lower bound to combined estimate (circle)"),
           plt.Line2D([], [], color=MUTED, lw=2, marker="s", ms=7, mec="white",
                      label="shuffled entrapment: lower bound to combined estimate (square)"),
           plt.Line2D([], [], color=INK, lw=0, marker="D", ms=6, mec="white",
                      label="shuffled entrapment: paired estimate (peptides, protein groups)"),
           plt.Line2D([], [], color=MUTED, lw=1, ls=(0, (4, 3)), label="nominal 1%")]
fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=7.5, frameon=False, labelcolor=MUTED)
fig.suptitle("Entrapment estimates of the false discovery proportion at a nominal 1% FDR (HYE, PXD028735)\n"
             "MaxQuant and FragPipe on Orbitrap DDA; DIA-NN on SWATH and on Orbitrap DIA of the same samples",
             fontsize=8.5, color=MUTED, x=0.01, ha="left")
fig.tight_layout(rect=(0, 0.08, 1, 0.94))
fig.savefig(HERE / "entrapment_fdr_hye_all.png", dpi=150)
print("wrote entrapment_fdr_hye_all.png")
