"""Figure for the CV post: median replicate CV by abundance quartile, per engine (from cv_by_quartile.csv)."""
import csv, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
rows = list(csv.DictReader(open("cv_by_quartile.csv")))
eng = ["MaxQuant raw Intensity", "FragPipe MaxLFQ", "MaxQuant LFQ", "DIA-NN pg_matrix"]
col = {"MaxQuant raw Intensity": "#b0b0b0", "MaxQuant LFQ": "#1f618d", "FragPipe MaxLFQ": "#c0392b", "DIA-NN pg_matrix": "#27ae60"}
fig, ax = plt.subplots(figsize=(8, 4.6), dpi=150)
for e in eng:
    r = [x for x in rows if x["engine"] == e]
    ys = [float(x["median_cv"]) for x in r]
    ax.plot(range(1, 5), ys, marker="o", lw=2.2, color=col[e], label=e)
    ax.annotate(f"{ys[0]:.1f}%", (1, ys[0]), textcoords="offset points", xytext=(-34, -4), fontsize=8, color=col[e])
ax.axhline(20, ls="--", lw=1, color="#888"); ax.text(4.05, 20.3, "20%", fontsize=8, color="#666")
ax.set_xticks(range(1, 5)); ax.set_xticklabels(["Q1\n(lowest)", "Q2", "Q3", "Q4\n(highest)"])
ax.set_xlabel("abundance quartile (mean intensity, condition A)"); ax.set_ylabel("median CV across 3 replicates (%)")
ax.set_title("Replicate CV falls with abundance, and depends on which column you use", fontsize=11, loc="left")
ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, fontsize=8); ax.set_ylim(0, 32)
fig.tight_layout(); fig.savefig(sys.argv[1] if len(sys.argv) > 1 else "cv_by_quartile.png")
