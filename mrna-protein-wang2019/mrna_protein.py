import numpy as np, pandas as pd
from scipy.stats import spearmanr

# Wang et al. 2019, Table EV6, sheet A: log10 FPKM (29 tissues), then log10 iBAQ (same 29 tissues)
t = pd.read_excel("Table_EV6.xlsx", sheet_name="A. Protein vs mRNA correlation", header=None)
mrna = t.iloc[2:, 2:31].astype(float).to_numpy()
prot = t.iloc[2:, 31:60].astype(float).to_numpy()

def rho(x, y):
    ok = ~np.isnan(x) & ~np.isnan(y)
    return spearmanr(x[ok], y[ok])[0] if ok.sum() >= 10 else np.nan

within_tissue = [rho(mrna[:, j], prot[:, j]) for j in range(29)]      # across genes, per tissue
within_gene = [rho(mrna[i], prot[i]) for i in range(len(mrna))]       # across tissues, per gene

print(f"{len(mrna)} genes x 29 tissues")
print(f"across genes, per tissue:  median rho {np.median(within_tissue):.2f} "
      f"(range {min(within_tissue):.2f}-{max(within_tissue):.2f})")
print(f"across tissues, per gene:  median rho {np.nanmedian(within_gene):.2f}, "
      f"{np.nanmean(np.array(within_gene) < 0):.0%} of genes negative")
