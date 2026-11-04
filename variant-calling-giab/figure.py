"""Figure for the HG002 chr20 benchmark, drawn from scores.csv (written by score.py).
Usage: python3 figure.py   -> variant_calling_giab.png next to this script (needs matplotlib)
  a, b  SNPs and indels: errors (false negatives + false positives) per call set and region, on a log scale
  c     recall vs precision of indels in all benchmark regions and in difficult regions
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
SETS = ["bcftools call -m (all calls)", "bcftools, QUAL >= 20", "GATK HaplotypeCaller (unfiltered)",
        "GATK HaplotypeCaller + hard filters", "DeepVariant"]
COLOUR = {"bcftools call -m (all calls)": "#86b6ef", "bcftools, QUAL >= 20": "#2a78d6",
          "GATK HaplotypeCaller (unfiltered)": "#f4a582", "GATK HaplotypeCaller + hard filters": "#eb6834",
          "DeepVariant": "#1baf7a"}
REGIONS = ["not in any difficult region", "tandem repeats and homopolymers", "segmental duplications", "low mappability",
           "all difficult regions", "all benchmark regions"]
SHORT = {"not in any difficult region": "easy", "tandem repeats and homopolymers": "TR / homopolymers",
         "segmental duplications": "segdups", "low mappability": "low mappability",
         "all difficult regions": "all difficult", "all benchmark regions": "all"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"

rows = list(csv.DictReader(open(HERE / "scores.csv")))
cov = max(int(r["dataset"].split("PCR-free ")[1].split("x")[0]) for r in rows if "PCR-free" in r["dataset"])
rows = [r for r in rows if f"PCR-free {cov}x" in r["dataset"]]
sets = [s for s in SETS if any(r["tool"] == s for r in rows)]


def val(s, vtype, region, metric):
    for r in rows:
        if r["tool"] == s and r["setting"] == f"{vtype}s" and r["subset"] == region and r["metric"] == metric:
            return float(r["value"])
    return None


fig, axs = plt.subplots(1, 3, figsize=(14, 4.4), gridspec_kw={"width_ratios": [5, 5, 3.4]})
w = 0.8 / len(sets)
for ax, vtype, title in ((axs[0], "SNP", "a  SNPs"), (axs[1], "INDEL", "b  Indels")):
    for k, s in enumerate(sets):
        ys = [(val(s, vtype, g, "fn") or 0) + (val(s, vtype, g, "fp") or 0) for g in REGIONS]
        ax.bar([i + (k - (len(sets) - 1) / 2) * w for i in range(len(REGIONS))], ys, w * 0.92, color=COLOUR[s], label=s)
    ax.set_yscale("log")
    ax.set_xticks(range(len(REGIONS)), [SHORT[g] for g in REGIONS], fontsize=7, rotation=20, ha="right")
    ax.set_title(f"{title}: errors (FN + FP) on chr20, HG002 {cov}x", fontsize=9, loc="left", color=INK)
    ax.set_ylabel("false negatives + false positives (log scale)", fontsize=8, color=MUTED)
    ax.tick_params(labelsize=7, colors=MUTED)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)
axs[0].legend(fontsize=6.5, frameon=False, loc="upper left")
ax = axs[2]
for s in sets:
    for region, marker in (("all benchmark regions", "o"), ("all difficult regions", "^")):
        r, p = val(s, "INDEL", region, "recall"), val(s, "INDEL", region, "precision")
        if r is not None and p is not None:
            ax.scatter(r, p, s=46, marker=marker, color=COLOUR[s], edgecolor="white", linewidth=0.8, zorder=3)
ax.scatter([], [], marker="o", color=MUTED, label="all benchmark regions")
ax.scatter([], [], marker="^", color=MUTED, label="all difficult regions")
ax.legend(fontsize=7, frameon=False, loc="lower left")
ax.set_xlabel("recall", fontsize=8, color=MUTED)
ax.set_ylabel("precision", fontsize=8, color=MUTED)
ax.set_title("c  Indels: recall vs precision", fontsize=9, loc="left", color=INK)
ax.tick_params(labelsize=7, colors=MUTED)
ax.grid(color=GRID, lw=0.6)
ax.set_axisbelow(True)
ax.spines[["top", "right"]].set_visible(False)
ax.spines[["left", "bottom"]].set_color(GRID)
fig.tight_layout()
fig.savefig(HERE / "variant_calling_giab.png", dpi=150)
print("wrote variant_calling_giab.png")
