"""Figure for the batch-correction study, drawn from scores.csv (written by batch_study.R).
Usage: python3 figure.py   -> batch_correction_unbalanced.png next to this script (needs matplotlib)
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
METHODS = ["ignore batch", "batch in the model", "ComBat (group protected)", "ComBat (batch only)",
           "removeBatchEffect (group protected)"]
STYLE = {"ignore batch": ("#999999", ":"), "batch in the model": ("#1b7837", "-"),
         "ComBat (group protected)": ("#b2182b", "-"), "ComBat (batch only)": ("#2166ac", "--"),
         "removeBatchEffect (group protected)": ("#ef8a62", "-.")}

rows = list(csv.DictReader(open(HERE / "scores.csv")))


def series(dataset_prefix, metric, subset=None):
    """{method: [(setting, value), ...]} in file order."""
    out = {m: [] for m in METHODS}
    for r in rows:
        if r["dataset"].startswith(dataset_prefix) and r["metric"] == metric and (subset is None or r["subset"] == subset):
            out[r["tool"]].append((r["setting"], float(r["value"])))
    return out


def k_label(setting):            # "batch 1 = 8 A + 4 B" -> "8:4"
    parts = setting.replace("batch 1 = ", "").replace(" A + ", ":").replace(" B", "")
    return parts


def panel(ax, data, label_of, title, ylabel, ref=None):
    for m in METHODS:
        xs = [label_of(s) for s, _ in data[m]]
        ys = [v for _, v in data[m]]
        c, ls = STYLE[m]
        ax.plot(range(len(xs)), ys, color=c, ls=ls, marker="o", ms=3.5, lw=1.6, label=m)
        ax.set_xticks(range(len(xs)), xs)
    if ref is not None:
        ax.axhline(ref, color="k", lw=0.7, ls="--")
    ax.set_title(title, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(labelsize=7.5)
    ax.spines[["top", "right"]].set_visible(False)


fig, axs = plt.subplots(2, 2, figsize=(9, 6.4))
panel(axs[0, 0], series("simulation (null)", "share_p_below_0.05"), k_label,
      "a  Simulation, no group effect", "share of genes with p < 0.05", 0.05)
panel(axs[0, 1], series("simulation (planted)", "power_bh_0.05"), k_label,
      "b  Simulation, 10% of genes changed: power", "changed genes found (BH 5%)")
panel(axs[1, 0], series("simulation (planted)", "fdp_bh_0.05"), k_label,
      "c  Simulation, 10% changed: false discovery proportion", "FDP of BH 5% calls", 0.05)
panel(axs[1, 1], series("bladderbatch", "share_p_below_0.05"), lambda s: s.split("= ")[1],
      "d  bladderbatch, fake labels within batches", "share of probes with p < 0.05", 0.05)
for ax in axs[:, 0].tolist() + [axs[0, 1]]:
    if ax is not axs[1, 1]:
        ax.set_xlabel("group A : group B in batch 1 (batch 2 reversed)", fontsize=8)
axs[1, 1].set_xlabel("share of fake group A in batches 1, 3, 5 (1 - that in 2, 4)", fontsize=8)
axs[0, 0].legend(fontsize=7, frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(HERE / "batch_correction_unbalanced.png", dpi=150)
print("wrote", HERE / "batch_correction_unbalanced.png")
