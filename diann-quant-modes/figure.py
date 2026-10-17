"""E. coli log2(B/A) by quantification route: MBR off (clean/nombr) and, with MBR on (clean/mbr), proteins that were
already scoreable without MBR versus proteins only scoreable with MBR. Same filters as score_modes.R / mbr_split.R."""
import re, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
S = ["A_Alpha", "A_Beta", "A_Gamma", "B_Alpha", "B_Beta", "B_Gamma"]
con = {l[1:].split()[0] for l in open("/home/sysoft/mq-hye/MaxQuant_v2.8.1.0/bin/conf/contaminants.fasta") if l.startswith(">")}
R = {"QuantUMS\n(default)": "pg_quantums_precision.tsv", "QuantUMS\nhigh-accuracy": "pg_quantums_accuracy.tsv",
     "Legacy\nMaxLFQ": "pg_legacy_direct.tsv", "Fragment\nMaxLFQ (iq)": "pg_fragment_maxlfq_iq.tsv",
     "Fragment\ndirectLFQ": "fragments_wide.aq_reformat.tsv.protein_intensities.tsv", "Precursor\nsum": "pg_precsum_quantums_precision.tsv",
     "Precursor sum,\nall-6-run precursors": "pg_precsum6_quantums_precision.tsv"}
def ecoli(d, f):
    try: t = pd.read_csv(f"clean/{d}/{f}", sep="\t")
    except FileNotFoundError: return None
    if "protein" in t.columns:
        nm = pd.read_csv(f"clean/{d}/pg_precsum_quantums_precision.tsv", sep="\t")[["Protein.Group", "Protein.Names"]]
        t = nm.merge(t.rename(columns={"protein": "Protein.Group"}), on="Protein.Group")
    t = t[~t["Protein.Group"].str.split(";").apply(lambda a: any(x in con for x in a))]
    sp = t["Protein.Names"].astype(str).apply(lambda x: set(re.findall(r"_(HUMAN|YEAST|ECOLI)\b", x)))
    t = t[sp.apply(lambda s: s == {"ECOLI"})]
    m = t[S].to_numpy(float); m[~np.isfinite(m) | (m <= 0)] = np.nan
    ok = (np.isfinite(m[:, :3]).sum(1) >= 2) & (np.isfinite(m[:, 3:]).sum(1) >= 2)
    lm = np.log2(m[ok]); return pd.Series(np.nanmean(lm[:, 3:], 1) - np.nanmean(lm[:, :3], 1), index=t["Protein.Group"][ok].to_numpy())
fig, ax = plt.subplots(figsize=(11, 4.6))
labels, data, cols = [], [], []
for i, (k, f) in enumerate(R.items()):
    off, on = ecoli("nombr", f), ecoli("mbr", f)
    if off is None: continue
    groups = [(off, "#4C72B0")]
    if on is not None:
        groups += [(on[on.index.isin(off.index)], "#55A868"), (on[~on.index.isin(off.index)], "#C44E52")]
    for j, (g, c) in enumerate(groups):
        pos = i * 4 + j
        b = ax.boxplot(g.to_numpy(), positions=[pos], widths=0.8, showfliers=False, patch_artist=True)
        b["boxes"][0].set_facecolor(c); b["medians"][0].set_color("black")
        ax.text(pos, -0.85, f"{100 * np.mean(np.abs(g - 2) <= 0.5):.0f}%", ha="center", fontsize=7)
    labels.append((i * 4 + 1, k))
ax.axhline(2, color="grey", ls="--", lw=1); ax.axhspan(1.5, 2.5, color="grey", alpha=0.08)
ax.set_xticks([p for p, _ in labels]); ax.set_xticklabels([k for _, k in labels], fontsize=8)
ax.set_ylim(-1.1, 5); ax.set_ylabel("E. coli log2(B/A)  (truth = 2)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#4C72B0", label="MBR off"), Patch(color="#55A868", label="MBR on: proteins also quantified without MBR"),
                   Patch(color="#C44E52", label="MBR on: proteins only quantified with MBR")], fontsize=7, loc="upper left")
ax.set_title("DIA-NN 2.3.2, PXD028735 SWATH, 3 vs 3: E. coli ratios by quantification route (outliers not drawn;\nnumbers = % within ±0.5 of truth)", fontsize=9)
plt.tight_layout(); plt.savefig("clean/ecoli_routes.png", dpi=150); print("saved")
