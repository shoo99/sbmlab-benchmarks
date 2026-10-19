"""Figure for the sbmlab post: string-hubs-measured.png"""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
h = json.load(open("hubs_results.json")); a = json.load(open("api_results.json"))
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
t = ["400", "700", "900"]
full = [h["pairs_by_threshold"][k] / 1000 for k in t]; notm = [h["pairs_by_threshold_no_textmining"][k] / 1000 for k in t]
x = np.arange(3); ax[0].bar(x - 0.2, full, 0.4, label="all channels", color="#3b6ea5"); ax[0].bar(x + 0.2, notm, 0.4, label="text mining removed", color="#e09f3e")
for i in range(3): ax[0].text(x[i] + 0.2, notm[i] + 15, f"-{100 - 100 * notm[i] / full[i]:.0f}%", ha="center", fontsize=9)
ax[0].set_xticks(x, ["≥0.4 (medium)", "≥0.7 (high)", "≥0.9 (highest)"]); ax[0].set_ylabel("human protein pairs (thousands)")
ax[0].set_title("A. Links kept at each cut-off"); ax[0].legend(frameon=False)
rng = np.random.default_rng(1)
for i, s in enumerate(["400", "700"]):
    v = np.array(a["summary"]["planted_chol_in_top10"][s]); ax[1].scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, color=["#888", "#2a9d8f"][i], s=28)
    ax[1].hlines(np.median(v), i - 0.25, i + 0.25, color="k")
ax[1].set_xticks([0, 1], ["cut-off 0.4", "cut-off 0.7"]); ax[1].set_ylim(-0.5, 10.5); ax[1].set_ylabel("cholesterol-synthesis proteins in top-10 hubs")
ax[1].set_title("B. 15 planted pathway proteins + 185 random:\nwho are the hubs?")
for i, k in enumerate(["s400_genome", "s400_bg"]):
    p = np.array([r[k]["p"] for r in a["random"]]); y = -np.log10(np.clip(p, 1e-30, 1))
    ax[2].scatter(i + rng.uniform(-0.12, 0.12, len(y)), y, color=["#c44536", "#3b6ea5"][i], s=28)
ax[2].axhline(-np.log10(0.05), ls="--", color="k", lw=1); ax[2].text(1.3, -np.log10(0.05) + 0.4, "p = 0.05", ha="right", fontsize=9)
ax[2].set_xticks([0, 1], ["whole-genome\nbackground (default)", "detected proteins\nas background"]); ax[2].set_ylabel("-log10 PPI enrichment p")
ax[2].set_title("C. 20 random lists of 200 detected proteins\n(cut-off 0.4)")
for g in ax: g.spines[["top", "right"]].set_visible(False)
fig.text(0.01, 0.01, "STRING 12.5, Homo sapiens. Detected proteins: DIA-NN on PXD028735 SWATH runs. sbmlab.com", fontsize=8, color="#555")
plt.tight_layout(rect=(0, 0.03, 1, 1)); plt.savefig("/home/sysoft/sbmlab-blog/public/images/blog/string-hubs-measured.png", dpi=130)
print("ok")
