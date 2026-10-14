import sys, pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
r = pd.read_csv("leak_k20.csv"); e = pd.read_csv("external_shuffled.csv")
s = r[r["labels"] == "shuffled"]
groups = [("No validation\n(training AUC)", s[s.method == "no CV (training AUC)"].auc, "#C44E52"),
          ("Leaky CV\n(selection before CV)", s[s.method == "leaky CV"].auc, "#C44E52"),
          ("Leaky 'external'\n(selection used test cohorts)", e.leaky, "#C44E52"),
          ("Proper CV\n(selection inside folds)", s[s.method == "proper CV"].auc, "#4C72B0"),
          ("Proper external\n(other cohorts untouched)", e.proper, "#4C72B0")]
fig, ax = plt.subplots(figsize=(10, 4.6))
b = ax.boxplot([g[1] for g in groups], widths=0.55, patch_artist=True, showfliers=False)
for p, g in zip(b["boxes"], groups): p.set_facecolor(g[2]); p.set_alpha(0.7)
rng = np.random.default_rng(0)
for i, g in enumerate(groups, 1): ax.scatter(i + rng.uniform(-0.15, 0.15, len(g[1])), g[1], s=8, color="k", alpha=0.5, zorder=3)
ax.axhline(0.5, ls="--", color="gray"); ax.text(0.55, 0.505, "true AUC = 0.5", ha="left", va="bottom", fontsize=8, color="gray")
ax.set_xticks(range(1, 6)); ax.set_xticklabels([g[0] for g in groups], fontsize=8.5)
ax.set_ylabel("AUC"); ax.set_ylim(0.3, 1.02)
ax.set_title("Shuffled labels (no real signal), 20 proteins chosen from 1,484, 50 permutations each", fontsize=10)
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=130); print("saved")
