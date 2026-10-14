import sys, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = pd.read_csv(sys.argv[1]); d["n"] = d[["Process", "Component", "Function", "KEGG"]].sum(axis=1)
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={"width_ratios": [4, 1.3]})
for ax, lst, title in ((axes[0], "random", "Random 300 detected proteins (20 draws)\ncorrect answer: nothing enriched"),
                       (axes[1], "planted", "15 cholesterol-synthesis\nproteins + 285 random (5 draws)")):
    s = d[d["list"] == lst]
    g = s[s.background == "genome"].set_index("seed")["n"]; t = s[s.background == "detected"].set_index("seed")["n"]
    x = range(len(g)); w = 0.4
    ax.bar([i - w/2 for i in x], g.values, w, label="background: whole genome (STRING default)", color="#C44E52")
    ax.bar([i + w/2 for i in x], t.values, w, label="background: detected proteins", color="#4C72B0")
    ax.set_xticks(list(x)); ax.set_xticklabels(g.index, fontsize=8); ax.set_xlabel("draw (seed)")
    ax.set_title(title, fontsize=10)
axes[0].set_ylabel("GO + KEGG terms at FDR < 0.05")
axes[0].legend(loc="upper left", fontsize=8, frameon=False)
fig.tight_layout(); fig.savefig(sys.argv[2], dpi=130); print("saved")
