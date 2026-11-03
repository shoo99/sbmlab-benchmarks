"""Figure for the microbiome DA study, drawn from scores.csv (written by da_study.R).
Usage: python3 figure.py   -> microbiome_da_methods.png next to this script (needs matplotlib)
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
rows = list(csv.DictReader(open(HERE / "scores.csv")))
METHODS = list(dict.fromkeys(r["tool"] for r in rows))


def value(scenario, metric, method):
    for r in rows:
        if r["dataset"].endswith(f"({scenario})") and r["metric"] == metric and r["tool"] == method:
            return float(r["value"])
    return float("nan")


fig, axs = plt.subplots(1, 3, figsize=(12, 4.2))
y = range(len(METHODS))
axs[0].barh(y, [value("null", "share_p_below_0.05", m) for m in METHODS], color="#777777")
axs[0].axvline(0.05, color="k", lw=0.8, ls="--")
axs[0].set_title("a  Nothing changed: share of OTUs with p < 0.05", fontsize=9)
if "ANCOM-BC2 + SS" in METHODS:
    i = METHODS.index("ANCOM-BC2 + SS")
    axs[0].text(0.02, i, "same p-values as ANCOM-BC2;\nthe filter only removes calls", va="center", fontsize=6.5, color="white")
for ax, metric, title, ref in [(axs[1], "fdp", "b  False discovery proportion of calls", 0.05),
                               (axs[2], "power", "c  Changed OTUs found (right direction)", None)]:
    h = 0.38
    ax.barh([i - h / 2 for i in y], [value("planted", metric, m) for m in METHODS], h, color="#4C72B0", label="15 OTUs changed 4-fold")
    ax.barh([i + h / 2 for i in y], [value("bloom", metric, m) for m in METHODS], h, color="#DD8452", label="same + top OTU 8-fold (bloom)")
    if ref is not None:
        ax.axvline(ref, color="k", lw=0.8, ls="--")
    if metric == "fdp":       # FDP from few calls is noisy: show the calls per repetition next to each bar
        for i, m in enumerate(METHODS):
            for dy, sc in ((-h / 2, "planted"), (h / 2, "bloom")):
                ax.text(value(sc, metric, m) + 0.01, i + dy, f"{value(sc, 'calls_per_rep', m):.1f} calls/rep", va="center", fontsize=6)
    ax.set_title(title, fontsize=9)
for ax in axs:
    ax.set_yticks(list(y), METHODS if ax is axs[0] else [""] * len(METHODS), fontsize=8)
    ax.invert_yaxis()
    ax.tick_params(axis="x", labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
axs[1].legend(fontsize=7, frameon=False, loc="lower right")
fig.suptitle("Microbiome differential abundance on smokers (left throat) 16S counts, 20 vs 20 samples, 20 repetitions per scenario", fontsize=9)
fig.tight_layout()
fig.savefig(HERE / "microbiome_da_methods.png", dpi=150)
print("wrote", HERE / "microbiome_da_methods.png")
