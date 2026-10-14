"""mRNA vs protein in 29 healthy human tissues (Wang et al. 2019, Mol Syst Biol 15:e8503, Table EV6 sheet A:
log10 FPKM and log10 iBAQ for genes detected as both in >= 10 tissues)."""
import sys, numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr

t = pd.read_excel(sys.argv[1], sheet_name="A. Protein vs mRNA correlation", header=None)
tissues = t.iloc[1, 2:31].tolist()
mrna = t.iloc[2:, 2:31].astype(float).to_numpy()
prot = t.iloc[2:, 31:60].astype(float).to_numpy()
genes = t.iloc[2:, 1].tolist()
paper_rho = t.iloc[2:, 60].astype(float).to_numpy()
print(f"genes: {len(genes)}, tissues: {len(tissues)}")

# 1. Within each tissue, across genes: does more mRNA mean more protein?
rows = []
for j, tis in enumerate(tissues):
    ok = ~np.isnan(mrna[:, j]) & ~np.isnan(prot[:, j])
    rho = spearmanr(mrna[ok, j], prot[ok, j])[0]
    r = pearsonr(mrna[ok, j], prot[ok, j])[0]
    ratio = prot[ok, j] - mrna[ok, j]                       # log10(protein/mRNA)
    q1, q3 = np.percentile(ratio, [25, 75]); lo, hi = np.percentile(ratio, [5, 95])
    rows.append(dict(tissue=tis, n=int(ok.sum()), spearman=rho, r2=r * r,
                     ratio_IQR_fold=10 ** (q3 - q1), ratio_90pct_fold=10 ** (hi - lo)))
w = pd.DataFrame(rows)
print("\nAcross genes, within each tissue:")
print(w[["n", "spearman", "r2", "ratio_IQR_fold", "ratio_90pct_fold"]].describe().loc[["min", "50%", "max"]].round(2))
print(w.sort_values("spearman")[["tissue", "n", "spearman", "r2"]].iloc[[0, 1, -2, -1]].round(2).to_string(index=False))

# 2. Within each gene, across tissues: does a tissue with more mRNA also have more protein?
rho_g = []
for i in range(len(genes)):
    ok = ~np.isnan(mrna[i]) & ~np.isnan(prot[i])
    rho_g.append(spearmanr(mrna[i, ok], prot[i, ok])[0] if ok.sum() >= 10 else np.nan)
rho_g = np.array(rho_g)
v = ~np.isnan(rho_g)
print(f"\nAcross tissues, within each gene ({v.sum()} genes): median Spearman {np.median(rho_g[v]):.2f}; "
      f"> 0.5: {np.mean(rho_g[v] > 0.5):.0%}; < 0: {np.mean(rho_g[v] < 0):.0%}")
print(f"agreement with the paper's own per-gene column: Pearson {pearsonr(rho_g[v], paper_rho[v])[0]:.3f}")
w.to_csv("within_tissue.csv", index=False); pd.DataFrame({"gene": genes, "rho_across_tissues": rho_g}).to_csv("per_gene.csv", index=False)

# 3. Significance of the per-gene correlations (paper's p-value column), variance effect, liver dynamic range
from scipy.stats import false_discovery_control
p_paper = t.iloc[2:, 61].astype(float).to_numpy()
print(f"negative & p<0.05: {int(((paper_rho < 0) & (p_paper < 0.05)).sum())}; negative & BH<0.05: "
      f"{int(((paper_rho < 0) & (false_discovery_control(p_paper) < 0.05)).sum())}; positive & p<0.05: {np.mean((paper_rho > 0) & (p_paper < 0.05)):.0%}")
q = pd.qcut(np.nanstd(mrna, axis=1), 4, labels=False)
print("median rho by quartile of mRNA SD across tissues:", [round(float(np.median(paper_rho[q == k])), 2) for k in range(4)])
j = tissues.index("Liver"); ok = ~np.isnan(mrna[:, j]) & ~np.isnan(prot[:, j])
span = lambda x: np.percentile(x, 99) - np.percentile(x, 1)
print(f"liver, {ok.sum()} genes with both: protein span {span(prot[ok, j]):.2f} log10, mRNA span {span(mrna[ok, j]):.2f} log10")
