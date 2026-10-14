"""Figures 1-5 for the HYE downstream-choices preprint. Inputs are the result tables in ../ (blog-pipeline/eval)."""
import gzip, io, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

E = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent / "fig"
OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "legend.fontsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 300})
SP = ["Human", "Yeast", "Ecoli"]; SPL = {"Human": "Human", "Yeast": "Yeast", "Ecoli": "E. coli"}
TRUTH = {"Human": 0, "Yeast": -1, "Ecoli": 2}
ENG = ["MaxQuant", "FragPipe", "DIA-NN"]; COL = {"MaxQuant": "#4C72B0", "FragPipe": "#DD8452", "DIA-NN": "#55A868"}


def label(ax, s):
    ax.text(-0.16, 1.13, s, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


# ---------- Figure 1: engines ----------
mq = pd.read_csv(E / "maxquant-hye/quant_scores.csv"); mq = mq[mq["quant"] == "LFQ intensity"].set_index("species")
fp = pd.read_csv(E / "fragpipe-hye/fragpipe_scores.csv"); fp = fp[fp["quant"] == "MaxLFQ Intensity"].set_index("species")
dn = pd.read_csv(E / "diann-hye/diann_scores.csv"); dn = dn[dn["set"].str.startswith("pg_matrix")].set_index("species")
common = pd.read_csv(E / "diann-hye/common_3866.csv")
tab = {"MaxQuant": mq, "FragPipe": fp, "DIA-NN": dn}

fig, axs = plt.subplots(2, 2, figsize=(7.2, 4.8)); axs = axs.ravel()
x = np.arange(3); w = 0.26
cen = pd.read_csv(Path(__file__).resolve().parent / "centred_engines.csv").set_index(["engine", "species"])
HATCH = {"DIA-NN": "///"}
for i, e in enumerate(ENG):
    lab = e + (" (SWATH, other instrument)" if e == "DIA-NN" else "")
    axs[0].bar(x + (i - 1) * w, [tab[e].loc[s, "n"] for s in SP], w, color=COL[e], label=lab, hatch=HATCH.get(e), edgecolor="white")
    axs[1].bar(x + (i - 1) * w, [100 * tab[e].loc[s, "within_0.5"] for s in SP], w, color=COL[e], hatch=HATCH.get(e), edgecolor="white")
    axs[1].scatter(x + (i - 1) * w, [cen.loc[(e, s), "within_centred"] for s in SP], s=12, color="k", zorder=3,
                   )
for a, t, yl in [(axs[0], "Proteins scored", "proteins"), (axs[1], "Within ±0.5 of truth, all scored\n(bars uncentred; dots centred on human median)", "% of proteins")]:
    a.set_xticks(x, [SPL[s] for s in SP]); a.set_title(t); a.set_ylabel(yl)
axs[1].set_ylim(0, 100); axs[0].set_ylim(0, 8600); axs[0].legend(frameon=False, loc="upper left", fontsize=6.5)
pos = []; data = []; cols = []
for j, s in enumerate(SP):
    sub = common[common["species"] == s]
    for i, (e, c) in enumerate(zip(ENG, ["fc_mq", "fc_fp", "fc_dn"])):
        pos.append(j * 4 + i); data.append(sub[c].values); cols.append(COL[e])
    axs[2].hlines(TRUTH[s], j * 4 - 0.5, j * 4 + 2.5, colors="k", linestyles="--", lw=0.7)
bp = axs[2].boxplot(data, positions=pos, widths=0.7, showfliers=False, patch_artist=True, medianprops={"color": "k", "lw": 0.8})
for p, c in zip(bp["boxes"], cols): p.set_facecolor(c); p.set_alpha(0.8)
axs[2].set_xticks([1, 5, 9], [f"{SPL[s]}\n(n={int((common['species'] == s).sum())})" for s in SP])
axs[2].set_ylabel("log2(B/A)"); axs[2].set_title("The 3,866 proteins all three quantified\n(dashed line = true ratio)")
mins = {"MaxQuant": 61.6, "FragPipe": 9.3 + 4.0, "DIA-NN": 377 + 11}
axs[3].barh(ENG, [mins[e] for e in ENG], color=[COL[e] for e in ENG], hatch=[HATCH.get(e, "") for e in ENG], edgecolor="white"); axs[3].set_xlim(0, 450)
axs[3].set_xlabel("wall time, minutes"); axs[3].set_title("Run time (6 runs)")
for e in ENG: axs[3].text(mins[e] + 6, e, f"{mins[e]:.0f} min", va="center", fontsize=7)
for a, l in zip(axs, "abcd"): label(a, l)
fig.tight_layout(h_pad=2.5); fig.savefig(OUT / "fig1_engines.png"); plt.close(fig)

# ---------- Figure 2: MBR ----------
pep = pd.read_csv(E / "maxquant-hye/mbr_scores.csv")
onoff = pd.read_csv(E / "maxquant-hye/mbr_onoff_LFQ_intensity.csv").set_index("species")
same = {"Human": (99.5, 99.1), "Yeast": (97.7, 95.0), "Ecoli": (94.8, 62.3)}  # maxquant-hye/README (compare_mbr_onoff.R)
dno = pd.read_csv(E / "diann-hye/diann_scores.csv")
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7))
grp = ["all 6 by MS/MS", "1-2 runs by matching", "3+ runs by matching"]; gl = ["all 6 MS/MS", "1–2 matched", "3+ matched"]
gc = ["#999999", "#8DA0CB", "#E78AC3"]
for i, g in enumerate(grp):
    v = [100 * pep[(pep.species == s) & (pep.group == g)]["within_0.5"].iloc[0] for s in SP]
    axs[0].bar(x + (i - 1) * w, v, w, color=gc[i], label=gl[i])
axs[0].set_ylim(0, 100); axs[0].set_xticks(x, [SPL[s] for s in SP]); axs[0].set_ylabel("% within ±0.5 (centred)")
axs[0].set_title("MaxQuant peptides\nby MBR origin"); axs[0].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=6)
axs[1].bar(x - w / 2, [same[s][0] for s in SP], w, color="#4C72B0", label="MBR on")
axs[1].bar(x + w / 2, [same[s][1] for s in SP], w, color="#BBBBBB", label="MBR off")
axs[1].set_ylim(0, 100); axs[1].set_xticks(x, [SPL[s] for s in SP]); axs[1].set_ylabel("% within ±0.5")
axs[1].set_title("MaxQuant LFQ proteins\nquantified with and without MBR"); axs[1].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=6)
dnon = dno[dno.set.str.startswith("pg_matrix")].set_index("species"); dnoff = dno[dno.set.str.startswith("no MBR")].set_index("species")
for i, (lab, on, off) in enumerate([("MaxQuant", onoff["n_on"], onoff["n_off"]), ("DIA-NN", dnon["n"], dnoff["n"])]):
    gain = [100 * (on[s] / off[s] - 1) for s in SP]
    axs[2].bar(x + (i - 0.5) * w, gain, w, color=COL[lab], label=lab)
axs[2].set_xticks(x, [SPL[s] for s in SP]); axs[2].set_ylabel("% more proteins scored with MBR")
axs[2].set_title("Proteins gained\nby MBR"); axs[2].set_ylabel("% more proteins scored"); axs[2].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=6)
for a, l in zip(axs, "abc"): label(a, l)
fig.tight_layout(); fig.savefig(OUT / "fig2_mbr.png"); plt.close(fig)

# ---------- Figure 3: imputation ----------
def parse(p):
    rows = []; data = filt = None
    for line in open(p):
        if line.startswith("===="):
            data = line.split("|")[0].replace("=", "").strip(); filt = "strict" if "BOTH" in line else "loose"
        elif "|" in line and not line.startswith("method"):
            f = [c.strip() for c in line.split("|")]
            found, of = f[5].split(" of ")
            rows.append(dict(data=data, filter=filt, method=f[0], ecoli_med=float(f[2]), found=int(found), of=int(of),
                             fp=int(f[6]), fdp=float(f[7].rstrip("%"))))
    return pd.DataFrame(rows)
imp = pd.concat([parse(E / "imputation-hye/impute_MBRon.out"), parse(E / "imputation-hye/impute_MBRoff.out")])
imp["recall"] = 100 * imp["found"] / imp["of"]
MC = {"none": "#000000", "MinDet": "#1b9e77", "MinProb": "#d95f02", "Perseus": "#7570b3", "QRILC": "#e7298a", "kNN": "#a6761d"}
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1, 1, 1.4]})
for a, d in zip(axs[:2], ["MBRon", "MBRoff"]):
    for _, r in imp[imp.data == d].iterrows():
        if r.method == "none" and r["filter"] == "loose": continue
        a.scatter(r.fdp, r.recall, color=MC[r.method], marker="o" if r["filter"] == "strict" else "^", s=22)
    a.set_xlabel("false discovery proportion (%)"); a.set_ylabel("true changes found (% of testable)")
    a.set_title(f"{'MBR on' if d == 'MBRon' else 'MBR off'}\n● strict filter  ▲ loose filter")
for m, c in MC.items(): axs[0].scatter([], [], color=c, label=m, s=18)
fig.legend(*axs[0].get_legend_handles_labels(), frameon=False, loc="lower center", ncol=6, bbox_to_anchor=(0.33, 0.0))
raw = pd.read_csv(io.BytesIO(gzip.open(E / "imputation-hye/impute_MBRoff.csv.gz").read()))
sub = raw[(raw.species == "Ecoli") & (raw["filter"].str.contains("EITHER"))]
meths = ["none", "MinDet", "MinProb", "Perseus", "QRILC", "kNN"]
vals = [sub[sub.method == m]["lfc_c"].dropna().values for m in meths]
bp = axs[2].boxplot(vals, showfliers=False, patch_artist=True, medianprops={"color": "k"})
for p, m in zip(bp["boxes"], meths): p.set_facecolor(MC[m]); p.set_alpha(0.6)
axs[2].axhline(2, ls="--", c="k", lw=0.7); axs[2].set_xticks(range(1, 7), meths, rotation=30)
axs[2].set_ylabel("E. coli log2(B/A), centred"); axs[2].set_title("MBR off, loose filter\nE. coli proteins (truth +2)")
for a, l in zip(axs, "abc"): label(a, l)
fig.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(OUT / "fig3_imputation.png"); plt.close(fig)

# ---------- Figure 4: statistics across engines ----------
st = pd.read_csv(Path(__file__).resolve().parent / "stats_three_engines.csv")
st["eng"] = st.engine.str.split().str[0]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.4), sharey=True)
tests = ["t", "limma", "DEqMS"]; xx = np.arange(3)
for a, norm in zip(axs, ["none", "median"]):
    for i, e in enumerate(ENG):
        s = st[(st.eng == e) & (st.norm == norm)]
        v0 = [s[(s.test == t) & (s.fc_min == 0)]["false_calls"].iloc[0] for t in tests]
        v5 = [s[(s.test == t) & (s.fc_min == 0.5)]["false_calls"].iloc[0] for t in tests]
        a.bar(xx + (i - 1) * w, v0, w, color=COL[e], alpha=0.85, label=e if norm == "none" else None)
        a.bar(xx + (i - 1) * w, v5, w * 0.5, color="k", alpha=0.6, label="with |log2FC| ≥ 0.5" if (norm == "none" and i == 0) else None)
    a.set_yscale("log"); a.set_xticks(xx, tests)
    a.set_title("Engine output as is" if norm == "none" else "After extra per-sample median centring")
axs[0].set_ylabel("unchanged human proteins called\n(adj. p < 0.05, log scale)")
h, l = axs[0].get_legend_handles_labels(); o = [l.index(k) for k in ENG + ["with |log2FC| ≥ 0.5"]]
axs[0].legend([h[i] for i in o], [l[i] for i in o], frameon=False, fontsize=6.5)
for a, l in zip(axs, "ab"): label(a, l)
fig.tight_layout(); fig.savefig(OUT / "fig4_statistics.png"); plt.close(fig)

# ---------- Figure 5: enrichment background ----------
cnt = pd.read_csv(E / "string-background/counts.csv")
cnt["total"] = cnt[["Process", "Component", "Function", "KEGG"]].sum(axis=1)
r = cnt[cnt.list == "random"]
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(6.4, 2.5))
g = r[r.background == "genome"]["total"].values; d = r[r.background == "detected"]["total"].values
for gi, di in zip(g, d): ax.plot([0, 1], [gi, di], c="#BBBBBB", lw=0.6)
ax.scatter(np.zeros_like(g), g, c="#C44E52", s=14); ax.scatter(np.ones_like(d), d, c="#4C72B0", s=14)
ax.set_xticks([0, 1], ["whole genome\n(STRING default)", "detected proteins"]); ax.set_xlim(-0.4, 1.4)
ax.set_ylabel("GO + KEGG terms at FDR < 0.05"); ax.set_title(f"20 random lists (no signal)\nmedian {np.median(g):.0f} vs {np.median(d):.0f}")
p = cnt[cnt.list == "planted"]; pg_ = p[p.background == "genome"]["total"].values; pd_ = p[p.background == "detected"]["total"].values
for gi, di in zip(pg_, pd_): ax2.plot([0, 1], [gi, di], c="#BBBBBB", lw=0.6)
ax2.scatter(np.zeros_like(pg_), pg_, c="#C44E52", s=14); ax2.scatter(np.ones_like(pd_), pd_, c="#4C72B0", s=14)
ax2.set_xticks([0, 1], ["whole genome", "detected proteins"]); ax2.set_xlim(-0.4, 1.4)
ax2.set_title("5 lists with 15 cholesterol-biosynthesis\nproteins planted (found in all 10 runs)")
label(ax, "a"); label(ax2, "b")
fig.tight_layout(); fig.savefig(OUT / "fig5_enrichment.png"); plt.close(fig)
print("ok", sorted(p.name for p in OUT.iterdir()))
