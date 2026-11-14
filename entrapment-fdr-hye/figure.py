"""Figure for the entrapment benchmark, drawn from scores.csv (written by score.py).
Usage: python3 figure.py   -> entrapment_fdr_hye.png next to this script (needs matplotlib)
One panel per level (precursors, peptides, protein groups, per run, MBR transfers only); one row per search.
The bar runs from the lower bound (fdp_lower) to the combined estimate (fdp_combined, the dot) at the nominal 1%
(DIA-NN q <= 0.01; MaxQuant and FragPipe as reported, MaxQuant protein groups at Q-value <= 0.01).
Filled dots: MBR on; open dots: MBR off. The dashed line is the nominal 1%.
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
