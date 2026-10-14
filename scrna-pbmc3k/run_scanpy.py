"""PBMC3k (10x, 2,700 cells) through the sbmlab single-cell guide's Scanpy steps, with real numbers.
Fixed QC thresholds from the guide (200-7500 genes, mt<20%) vs per-dataset MAD thresholds."""
import time, resource, sys, numpy as np, scanpy as sc, matplotlib
matplotlib.use("Agg")
sc.settings.verbosity = 0; sc.settings.figdir = sys.argv[2]
T = {}; t0 = time.time()
adata = sc.read_10x_mtx(sys.argv[1], var_names="gene_symbols", cache=False)
adata.var_names_make_unique()
print("raw:", adata.n_obs, "cells x", adata.n_vars, "genes")
sc.pp.filter_cells(adata, min_genes=200); sc.pp.filter_genes(adata, min_cells=3)
print("after min_genes=200 / min_cells=3:", adata.n_obs, "x", adata.n_vars)
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True)
q = adata.obs
print(f"genes/cell median {q.n_genes_by_counts.median():.0f} (max {q.n_genes_by_counts.max():.0f}); mt% median {q.pct_counts_mt.median():.2f} (99th pct {q.pct_counts_mt.quantile(.99):.1f})")
fixed = (q.n_genes_by_counts < 7500) & (q.pct_counts_mt < 20)
print("fixed thresholds (genes<7500, mt<20%) keep:", int(fixed.sum()), "of", adata.n_obs)
def mad_out(x, n):
    med = np.median(x); mad = np.median(np.abs(x - med)); return (x < med - n * mad) | (x > med + n * mad)
mad = ~(mad_out(q.log1p_total_counts, 5) | mad_out(q.log1p_n_genes_by_counts, 5) | (q.pct_counts_mt > np.median(q.pct_counts_mt) + 3 * np.median(np.abs(q.pct_counts_mt - np.median(q.pct_counts_mt)))))
print("MAD thresholds (5 MAD counts/genes, 3 MAD mt%) keep:", int(mad.sum()), "of", adata.n_obs, "| mt% cut at", round(float(np.median(q.pct_counts_mt) + 3 * np.median(np.abs(q.pct_counts_mt - np.median(q.pct_counts_mt)))), 2))
seu = (q.n_genes_by_counts > 200) & (q.n_genes_by_counts < 2500) & (q.pct_counts_mt < 5)
print("Seurat PBMC3k tutorial thresholds (200<genes<2500, mt<5%) keep:", int(seu.sum()))
print("cells removed by MAD but kept by fixed:", int((fixed & ~mad).sum()), "| removed by fixed but kept by MAD:", int((~fixed & mad).sum()))
sc.pp.scrublet(adata)
print(f"scrublet predicted doublets: {int(adata.obs.predicted_doublet.sum())} of {adata.n_obs} ({100*adata.obs.predicted_doublet.mean():.1f}%)")
adata = adata[mad & ~adata.obs.predicted_doublet].copy()
print("after MAD QC + doublet removal:", adata.n_obs, "cells"); T["load+QC"] = time.time() - t0; t0 = time.time()
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", layer="counts")
sc.tl.pca(adata, n_comps=50, mask_var="highly_variable")
vr = adata.uns["pca"]["variance_ratio"]; print("PCA variance explained PC1-10:", np.round(vr[:10]*100, 1).tolist())
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30); sc.tl.umap(adata)
for res in (0.5, 1.0):
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}", flavor="igraph", n_iterations=2, random_state=0)
    print(f"leiden resolution {res}: {adata.obs[f'leiden_{res}'].nunique()} clusters")
T["norm..cluster"] = time.time() - t0; t0 = time.time()
markers = {"T cells": ["CD3E","IL7R"], "CD8 T": ["CD8A"], "NK": ["GNLY","NKG7"], "B cells": ["MS4A1","CD79A"],
           "CD14 mono": ["CD14","LYZ"], "FCGR3A mono": ["FCGR3A","MS4A7"], "DC": ["FCER1A","CST3"], "Platelets": ["PPBP"]}
sc.tl.rank_genes_groups(adata, "leiden_0.5", method="wilcoxon")
top = {c: [adata.uns["rank_genes_groups"]["names"][c][i] for i in range(5)] for c in adata.obs["leiden_0.5"].cat.categories}
for c, g in top.items(): print(f"cluster {c} (n={int((adata.obs['leiden_0.5']==c).sum())}): {', '.join(g)}")
T["markers"] = time.time() - t0
sc.pl.umap(adata, color=["leiden_0.5", "CD3E", "MS4A1", "LYZ", "FCGR3A", "GNLY"], ncols=3, save="_pbmc3k.png", show=False)
sc.pl.dotplot(adata, markers, groupby="leiden_0.5", save="_pbmc3k_markers.png", show=False)
print("timings (s):", {k: round(v, 1) for k, v in T.items()}, "| peak RSS GB:", round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 2))
adata.write_h5ad(sys.argv[2] + "/pbmc3k_processed.h5ad")
