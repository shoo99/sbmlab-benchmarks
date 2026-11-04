"""Figure for the mock-community benchmark, drawn from scores.csv (written by score.py).
Usage: python3 figure.py   -> metagenome_profilers_mock.png next to this script (needs matplotlib)
  a  D6300 (even): estimated vs true percentage of each mock species, per profile
  b  D6310 (log distribution): estimated vs true percentage on log scales; missing species drawn at the bottom
  c  L1 distance to the true composition per profile and data set
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
rows = list(csv.DictReader(open(HERE / "scores.csv")))
truth = {}
for r in csv.DictReader(open(HERE / "truth.tsv"), delimiter="\t"):
    truth.setdefault(r["dataset"], {})[r["species"]] = float(r["gdna_percent"])
tools = []
for r in rows:
    if r["metric"] == "estimated_percent" and r["tool"] not in tools:
        tools.append(r["tool"])
colour = {t: PALETTE[i % len(PALETTE)] for i, t in enumerate(tools)}


def est(ds_start, tool, species):
    for r in rows:
        if r["dataset"].startswith(ds_start) and r["tool"] == tool and r["subset"] == species and r["metric"] == "estimated_percent":
            return float(r["value"])
    return None


def style(ax, title, ylabel):
    ax.set_title(title, fontsize=9, loc="left", color=INK)
    ax.set_ylabel(ylabel, fontsize=8, color=MUTED)
    ax.tick_params(labelsize=7, colors=MUTED)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)


fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), gridspec_kw={"width_ratios": [6, 4, 3]})
ax = axs[0]
sp = list(truth["D6300"])
tl = [t for t in tools if est("ZymoBIOMICS D6300", t, sp[0]) is not None]
w = 0.8 / max(len(tl), 1)
for k, t in enumerate(tl):
    ax.bar([i + (k - (len(tl) - 1) / 2) * w for i in range(len(sp))], [est("ZymoBIOMICS D6300", t, s) for s in sp],
           w * 0.92, color=colour[t], label=t)
for i, s in enumerate(sp):
    ax.plot([i - 0.45, i + 0.45], [truth["D6300"][s]] * 2, color=INK, lw=1.2)
ax.set_xticks(range(len(sp)), [s.split()[0][0] + ". " + s.split()[1] for s in sp], fontsize=7, rotation=25, ha="right")
style(ax, "a  D6300 even: estimated % per species (black line: true genomic-DNA %)", "estimated % of reads")
ax.legend(fontsize=7, frameon=False, loc="upper right")
ax = axs[1]
floor = 1e-5
shift = {t: 10 ** (0.08 * (k - (len(tools) - 1) / 2)) for k, t in enumerate(tools)}   # small horizontal offsets
markers = ["o", "D", "s", "^", "v", "P"]
for k, t in enumerate(tools):
    xs, ys = [], []
    for s, tv in truth["D6310"].items():
        e = est("ZymoBIOMICS D6310", t, s)
        if e is None:
            continue
        xs.append(tv * shift[t])
        ys.append(max(e, floor))
    if xs:
        ax.scatter(xs, ys, s=34, marker=markers[k % len(markers)], color=colour[t], edgecolor="white", linewidth=0.8,
                   zorder=3, label=t)
lims = [floor, 200]
ax.plot(lims, lims, color=INK, lw=0.8, ls="--")
ax.axhline(floor, color=GRID, lw=0.8)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(5e-5, 200)
ax.set_ylim(floor / 2, 200)
ax.set_xlabel("true genomic-DNA %", fontsize=8, color=MUTED)
style(ax, "b  D6310 log: estimated vs true (bottom line: not detected)", "estimated % (log)")
ax.legend(fontsize=6.5, frameon=False, loc="upper left")
ax = axs[2]
sets = [("ZymoBIOMICS D6300", "D6300"), ("ZymoBIOMICS D6310", "D6310")]
w = 0.8 / max(len(tools), 1)
for k, t in enumerate(tools):
    ys = []
    for start, _ in sets:
        v = [float(r["value"]) for r in rows if r["dataset"].startswith(start) and r["tool"] == t and r["metric"] == "l1_distance"]
        ys.append(v[0] if v else 0)
    ax.bar([i + (k - (len(tools) - 1) / 2) * w for i in range(len(sets))], ys, w * 0.92, color=colour[t])
ax.set_xticks(range(len(sets)), [s for _, s in sets], fontsize=8)
style(ax, "c  L1 distance to the truth", "sum of |estimated - true| (%)")
fig.tight_layout()
fig.savefig(HERE / "metagenome_profilers_mock.png", dpi=150)
print("wrote metagenome_profilers_mock.png")
