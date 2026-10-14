import sys, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = pd.read_csv(sys.argv[1]); out = sys.argv[2]
tools = [("fc_mq", "MaxQuant\n(Orbitrap DDA)"), ("fc_fp", "FragPipe\n(Orbitrap DDA)"), ("fc_dn", "DIA-NN\n(TripleTOF SWATH)")]
sp = [("Human", 0, "#4C72B0"), ("Yeast", -1, "#DD8452"), ("Ecoli", 2, "#55A868")]
fig, ax = plt.subplots(figsize=(9, 5))
pos, data, cols = [], [], []
for i, (c, lab) in enumerate(tools):
    for j, (s, t, col) in enumerate(sp):
        pos.append(i * 4 + j); data.append(d.loc[d.species == s, c].values); cols.append(col)
b = ax.boxplot(data, positions=pos, widths=0.8, patch_artist=True, showfliers=False)
for p, col in zip(b["boxes"], cols): p.set_facecolor(col); p.set_alpha(0.75)
for s, t, col in sp: ax.axhline(t, ls="--", lw=1, color=col)
ax.set_xticks([1, 5, 9]); ax.set_xticklabels([l for _, l in tools])
ax.set_ylabel("log2(B/A)"); ax.set_title("Same 3,866 proteins quantified by all three tools (dashed = true ratio)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=c, alpha=0.75, label=f"{s} (truth {t:+d})") for s, t, c in sp], loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)
fig.tight_layout(); fig.savefig(out, dpi=130); print("saved", out)
