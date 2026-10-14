import sys, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
E = sys.argv[1]; out = sys.argv[2]
d = pd.concat([pd.read_csv(f"{E}/impute_MBRon.csv.gz"), pd.read_csv(f"{E}/impute_MBRoff.csv.gz")])
d = d[d["filter"] == ">=2 of 3 in EITHER group"]
order = ["none", "MinDet", "MinProb", "Perseus", "QRILC", "kNN"]
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), dpi=150, sharey=True)
for i, tag in enumerate(["MBRoff", "MBRon"]):
    x = d[(d.data == tag) & (d.species == "Ecoli")]
    data = [x[x.method == m].lfc_c.values for m in order]
    ax[i].boxplot(data, tick_labels=order, showfliers=False)
    ax[i].axhline(2, color="#d1495b", ls="--", lw=1); ax[i].text(6.4, 2.05, "truth +2", color="#d1495b", fontsize=8, ha="right")
    ax[i].axhline(0, color="grey", lw=0.5)
    n = [len(v) for v in data]
    ax[i].set_title(f"E. coli proteins, MaxQuant {'without' if tag=='MBRoff' else 'with'} MBR (n={n[1]}; 'none' n={n[0]})", fontsize=9)
    ax[i].tick_params(axis="x", labelsize=8.5)
ax[0].set_ylabel("log2(B/A), centred on human median")
fig.suptitle("Filter '≥2 of 3 in either group': how each imputation method recovers a known +2 change", fontsize=10)
fig.tight_layout(); fig.savefig(out); print("ok")
