"""GLUE (Cao & Gao 2022, Nat Biotechnol) on the unpaired RNA and ATAC of the 10x Multiome PBMC data (issue #16).

Follows the scglue tutorial: RNA counts with 2,000 highly variable genes (Seurat v3 flavour), normalised, log, scaled,
100 PCs; ATAC LSI (100 components); an RNA-anchored guidance graph from gene coordinates (gene body plus promoter,
from work/gene_coords.tsv written by seurat_methods.R) to peaks; a SCGLUE model with negative-binomial likelihoods,
trained on CPU (random seed 0). The ATAC cells get new names, so the model cannot pair them by barcode.
The cells and their RNA cell types are those of seurat_methods.R (work/cells.tsv).
Scoring as in seurat_methods.R: FOSCTTM in both directions in the GLUE embedding, and labels from the majority of the
15 nearest RNA cells. Writes results/glue.tsv and prints a summary.
Usage: python run_glue.py [data_dir]     (needs scglue, torch, scanpy, h5py; 1-3 h on 4 CPU threads; a rerun resumes
       from the last saved training stage in work/glue)
"""
import sys
import time
from itertools import chain
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scglue
import torch

HERE = Path(__file__).resolve().parent
data = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "data"
work, res = HERE / "work", HERE / "results"
torch.set_num_threads(4)
np.random.seed(0)
t0 = time.time()

with h5py.File(data / "pbmc10k_multiome.h5", "r") as f:
    m = f["matrix"]
    X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=m["shape"][:])
    ftype = m["features/feature_type"][:].astype(str)
    fname = m["features/name"][:].astype(str)
    fint = m["features/interval"][:].astype(str)
    barcodes = m["barcodes"][:].astype(str)
cells = pd.read_csv(work / "cells.tsv", sep="\t")
col = pd.Series(np.arange(len(barcodes)), index=barcodes)[cells["barcode"]].to_numpy()
X = X[:, col]
g = ftype == "Gene Expression"
p = ftype == "Peaks"
rna = ad.AnnData(sp.csr_matrix(X[g].T, dtype=np.float32), obs=pd.DataFrame(index=cells["barcode"].to_numpy()),
                 var=pd.DataFrame(index=pd.Index(fname[g]).astype(str)))
rna.var_names_make_unique()
atac = ad.AnnData(sp.csr_matrix(X[p].T, dtype=np.float32),
                  obs=pd.DataFrame(index=[f"atac-{i}" for i in range(X.shape[1])]),     # no barcode to pair by
                  var=pd.DataFrame(index=fint[p]))
del X
coords = pd.read_csv(work / "gene_coords.tsv", sep="\t").drop_duplicates("gene").set_index("gene")
rna = rna[:, rna.var_names.isin(coords.index)].copy()
rna.var["chrom"] = coords.loc[rna.var_names, "chr"].to_numpy()
rna.var["chromStart"] = coords.loc[rna.var_names, "start"].to_numpy() - 1
rna.var["chromEnd"] = coords.loc[rna.var_names, "end"].to_numpy()
rna.var["strand"] = coords.loc[rna.var_names, "strand"].to_numpy()
split = atac.var_names.str.split(r"[:-]", regex=True)
atac.var["chrom"] = split.str[0]
atac.var["chromStart"] = split.str[1].astype(int)
atac.var["chromEnd"] = split.str[2].astype(int)
print(f"== {rna.n_obs} cells; RNA {rna.n_vars} genes with coordinates; ATAC {atac.n_vars} peaks", flush=True)

rna.layers["counts"] = rna.X.copy()
sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto", random_state=0)
scglue.data.lsi(atac, n_components=100, n_iter=15, random_state=0)
guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
scglue.graph.check_graph(guidance, [rna, atac])
scglue.models.configure_dataset(rna, "NB", use_highly_variable=True, use_layer="counts", use_rep="X_pca")
scglue.models.configure_dataset(atac, "NB", use_highly_variable=True, use_rep="X_lsi")
hv = guidance.subgraph(chain(rna.var.query("highly_variable").index, atac.var.query("highly_variable").index)).copy()
print(f"   guidance graph: {guidance.number_of_nodes()} nodes, {guidance.number_of_edges()} edges; "
      f"{int(atac.var['highly_variable'].sum())} peaks linked to the variable genes; setup {time.time() - t0:.0f} s", flush=True)
# The two training stages of scglue.models.fit_SCGLUE (pretraining, then fine-tuning with the discriminator and
# balancing weights), written out so that each stage is saved and a rerun resumes from the last saved stage.
# Each stage can take up to 152 epochs of about 35 s on 4 CPU threads.
adatas = {"rna": rna, "atac": atac}
gdir = work / "glue"
pre_file, fine_file = gdir / "pretrain" / "pretrain.dill", gdir / "fine-tune" / "fine-tune.dill"
if fine_file.exists():
    glue = scglue.models.load_model(fine_file)
    print(f"   loaded the fine-tuned model from {fine_file.relative_to(HERE)}", flush=True)
else:
    if pre_file.exists():
        pretrain = scglue.models.load_model(pre_file)
        print(f"   loaded the pretrained model from {pre_file.relative_to(HERE)}", flush=True)
    else:
        pretrain = scglue.models.SCGLUEModel(adatas, sorted(hv.nodes), random_seed=0, shared_batches=False)
        pretrain.compile()
        pretrain.fit(adatas, hv, align_burnin=np.inf, safe_burnin=False, directory=str(gdir / "pretrain"))
        pretrain.save(pre_file)
        print(f"   pretrained in {time.time() - t0:.0f} s", flush=True)
    for k, a in adatas.items():
        a.obsm["X_pretrain"] = pretrain.encode_data(k, a)
    scglue.data.estimate_balancing_weight(rna, atac, use_rep="X_pretrain", key_added="balancing_weight")
    for a in adatas.values():
        scglue.data.get_dataset_config(a)["use_dsc_weight"] = "balancing_weight"
        del a.obsm["X_pretrain"]
    glue = scglue.models.SCGLUEModel(adatas, sorted(hv.nodes), random_seed=0)
    glue.adopt_pretrained_model(pretrain)
    glue.compile()
    glue.random_seed += 1                      # as fit_SCGLUE: a different data order from pretraining
    glue.fit(adatas, hv, directory=str(gdir / "fine-tune"))
    glue.save(fine_file)
R = glue.encode_data("rna", rna)
A = glue.encode_data("atac", atac)
print(f"   trained in {time.time() - t0:.0f} s; embedding {R.shape[1]} dimensions", flush=True)


def foscttm(A, R):
    """Per-cell FOSCTTM both ways; row i of A and of R are the same nucleus."""
    na, nr = (A ** 2).sum(1), (R ** 2).sum(1)
    dtrue = na + nr - 2 * (A * R).sum(1)
    n = len(A)
    fa, fr = np.empty(n), np.empty(n)
    for s in range(0, n, 1000):
        e = min(n, s + 1000)
        d = na[s:e, None] + nr[None, :] - 2 * A[s:e] @ R.T
        fa[s:e] = (d < dtrue[s:e, None]).sum(1) / (n - 1)
        d2 = nr[s:e, None] + na[None, :] - 2 * R[s:e] @ A.T
        fr[s:e] = (d2 < dtrue[s:e, None]).sum(1) / (n - 1)
    return fa, fr


def knn_labels(A, R, labels, k=15):
    nr = (R ** 2).sum(1)
    out = []
    for s in range(0, len(A), 1000):
        d = (A[s:s + 1000] ** 2).sum(1)[:, None] + nr[None, :] - 2 * A[s:s + 1000] @ R.T
        for row in np.argsort(d, axis=1)[:, :k]:
            v, c = np.unique(labels[row], return_counts=True)
            out.append(v[np.argmax(c)])
    return np.array(out)


A64, R64 = A.astype(np.float64), R.astype(np.float64)
fa, fr = foscttm(A64, R64)
true = cells["cell_type"].to_numpy()
pred = knn_labels(A64, R64, true)
out = pd.DataFrame({"barcode": cells["barcode"], "true_type": true, "predicted_type": pred,
                    "foscttm_atac": np.round(fa, 6), "foscttm_rna": np.round(fr, 6)})
out.to_csv(res / "glue.tsv", sep="\t", index=False)
print(f"   glue: FOSCTTM {np.mean(np.r_[fa, fr]):.4f} (ATAC->RNA {fa.mean():.4f}, RNA->ATAC {fr.mean():.4f}), "
      f"label accuracy {np.mean(pred == true):.3f}; total {time.time() - t0:.0f} s")
from importlib.metadata import version  # noqa: E402
print(f"Python {sys.version.split()[0]} | scglue {version('scglue')} | torch {torch.__version__} | "
      f"scanpy {version('scanpy')} | numpy {np.__version__}")
