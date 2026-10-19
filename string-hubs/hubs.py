"""STRING 12.5 human network: what the confidence score is made of, and whether hubs are biology or attention.

Inputs (downloaded 2026-10-19 into /home/sysoft/data/string125, not in git):
  9606.protein.links.detailed.v12.5.txt.gz, 9606.protein.info.v12.5.txt.gz (stringdb-downloads.org)
  gene2pubmed.gz, Homo_sapiens.gene_info.gz (NCBI)
  CEGv2.txt, NEGv1.txt (Hart et al. 2017, hart-lab/bagel)
  rna_celline.tsv.zip (Human Protein Atlas cell line RNA, file dated 2025-11-05)
Combined score recomputed with STRING's published procedure (FAQ / combine_subscores.v2.py): remove prior 0.041 from
each channel, combine as 1 - prod(1 - s), add the prior back once. 'No text mining' = same with the textmining channel set to 0.
Outputs: hubs_results.json, degree_table.tsv.gz
"""
import json, gzip, sys
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm

D = sys.argv[1] if len(sys.argv) > 1 else "/home/sysoft/data/string125"
CH = ["neighborhood", "fusion", "cooccurence", "coexpression", "experimental", "database", "textmining"]
PRIOR = 0.041
out = {}

L = pd.read_csv(f"{D}/9606.protein.links.detailed.v12.5.txt.gz", sep=" ",
                dtype={c: "int16" for c in CH + ["combined_score"]})
L = L[L.protein1 < L.protein2].reset_index(drop=True)          # each pair is listed twice
out["n_pairs_total"] = int(len(L))

def combine(df, chans):
    one_minus = np.ones(len(df))
    for c in chans:
        s = df[c].to_numpy() / 1000.0
        s = np.where(s < PRIOR, PRIOR, s)
        one_minus *= 1.0 - (s - PRIOR) / (1 - PRIOR)
    tot = (1.0 - one_minus) * (1 - PRIOR) + PRIOR
    return np.floor(tot * 1000).astype(int)

L["recomp"] = combine(L, CH)
diff = (L.recomp - L.combined_score).abs()
out["recompute_check"] = {"exact": round(100 * float((diff == 0).mean()), 2), "within_1": round(100 * float((diff <= 1).mean()), 2),
                          "within_5": round(100 * float((diff <= 5).mean()), 2), "max_abs_diff": int(diff.max())}
L["no_tm"] = combine(L, [c for c in CH if c != "textmining"])

THR = [150, 400, 700, 900]
out["pairs_by_threshold"] = {t: int((L.combined_score >= t).sum()) for t in THR}
out["pairs_by_threshold_no_textmining"] = {t: int((L.no_tm >= t).sum()) for t in THR}

# which channels carry each edge at a threshold
def channel_profile(t):
    e = L[L.combined_score >= t]
    nz = e[CH] > 0
    only_tm = nz.textmining & ~nz[[c for c in CH if c != "textmining"]].any(axis=1)
    falls = e.no_tm < t                                           # edge drops below threshold without text mining
    return {"n": int(len(e)),
            "has_channel_pct": {c: round(100 * float(nz[c].mean()), 1) for c in CH},
            "textmining_only_pct": round(100 * float(only_tm.mean()), 1),
            "below_threshold_without_textmining_pct": round(100 * float(falls.mean()), 1),
            "has_experimental_or_database_pct": round(100 * float((nz.experimental | nz.database).mean()), 1)}
out["channels"] = {t: channel_profile(t) for t in (400, 700, 900)}

# degree
def degree(col, t):
    e = L[L[col] >= t]
    return pd.concat([e.protein1, e.protein2]).value_counts()

info = pd.read_csv(f"{D}/9606.protein.info.v12.5.txt.gz", sep="\t").rename(columns={"#string_protein_id": "sid"})
deg = pd.DataFrame({"sid": info.sid, "symbol": info.preferred_name})
for col, t, name in (("combined_score", 400, "deg400"), ("combined_score", 700, "deg700"),
                     ("no_tm", 400, "deg400_notm"), ("no_tm", 700, "deg700_notm")):
    deg[name] = deg.sid.map(degree(col, t)).fillna(0).astype(int)

# symbol -> GeneID (official symbol first, then unique synonym), GeneID -> number of PubMed records
gi = pd.read_csv(f"{D}/Homo_sapiens.gene_info.gz", sep="\t", usecols=["#tax_id", "GeneID", "Symbol", "Synonyms", "type_of_gene"])
sym2id = dict(zip(gi.Symbol, gi.GeneID))
syn = gi.assign(s=gi.Synonyms.str.split("|")).explode("s")
syn = syn[syn.s.notna() & (syn.s != "-")]
syn_unique = syn.groupby("s").GeneID.nunique()
syn2id = dict(syn[syn.s.isin(syn_unique[syn_unique == 1].index)][["s", "GeneID"]].drop_duplicates().values)
deg["GeneID"] = deg.symbol.map(sym2id)
miss = deg.GeneID.isna()
deg.loc[miss, "GeneID"] = deg.loc[miss, "symbol"].map(syn2id)
out["mapping"] = {"proteins": int(len(deg)), "mapped_to_geneid": int(deg.GeneID.notna().sum())}

pubs = {}
with gzip.open(f"{D}/gene2pubmed.gz", "rt") as f:
    next(f)
    for line in f:
        if line.startswith("9606\t"):
            g = int(line.split("\t", 2)[1]); pubs[g] = pubs.get(g, 0) + 1
deg["pubs"] = deg.GeneID.map(lambda g: pubs.get(int(g), 0) if pd.notna(g) else np.nan)
m = deg[deg.GeneID.notna()].copy()

def spearman(a, b):
    r = stats.spearmanr(a, b); return round(float(r.statistic), 3)
out["spearman_degree_pubs"] = {k: spearman(m[k], m.pubs) for k in ("deg400", "deg700", "deg400_notm", "deg700_notm")}

# top hubs and how well-studied they are
m["pubs_pct"] = m.pubs.rank(pct=True)
def hubs(col, k=100):
    h = m.nlargest(k, col)
    return {"median_pubs": float(h.pubs.median()), "median_pubs_pct": round(100 * float(h.pubs_pct.median()), 1),
            "top10": h.head(10)[["symbol", col, "pubs"]].values.tolist()}
out["top_hubs"] = {c: hubs(c) for c in ("deg400", "deg700", "deg700_notm")}
out["all_proteins_median_pubs"] = float(m.pubs.median())

# essential (CEGv2) vs non-essential (NEGv1): degree alone, publications alone, both
ceg = pd.read_csv(f"{D}/CEGv2.txt", sep="\t"); neg = pd.read_csv(f"{D}/NEGv1.txt", sep="\t")
lab = pd.concat([pd.DataFrame({"GeneID": ceg.ENTREZ_ID, "ess": 1}), pd.DataFrame({"GeneID": neg.ENTREZ_ID, "ess": 0})])
e = m.merge(lab, on="GeneID").drop_duplicates("GeneID")
out["ess_n"] = {"essential": int(e.ess.sum()), "nonessential": int((e.ess == 0).sum())}
def auc(score, y):
    r = stats.rankdata(score); n1 = y.sum(); n0 = len(y) - n1
    return round(float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)), 3)
y = e.ess.to_numpy()
out["auc_essential"] = {c: auc(e[c].to_numpy(), y) for c in ("pubs", "deg400", "deg700", "deg400_notm", "deg700_notm")}
out["median_by_ess"] = {c: {"essential": float(e[e.ess == 1][c].median()), "nonessential": float(e[e.ess == 0][c].median())}
                        for c in ("pubs", "deg700", "deg700_notm")}
# logistic: essential ~ log1p(degree) + log1p(pubs), standardised
lr = {}
for c in ("deg700", "deg700_notm"):
    X = pd.DataFrame({"degree": np.log1p(e[c]), "pubs": np.log1p(e.pubs)})
    X = (X - X.mean()) / X.std()
    fit = sm.Logit(y, sm.add_constant(X)).fit(disp=0)
    lr[c] = {"coef_degree": round(float(fit.params["degree"]), 3), "p_degree": float(fit.pvalues["degree"]),
             "coef_pubs": round(float(fit.params["pubs"]), 3), "p_pubs": float(fit.pvalues["pubs"])}
out["logit_essential"] = lr
# within publication tertiles: is degree still higher for essential genes?
e["pub_bin"] = pd.qcut(e.pubs.rank(method="first"), 3, labels=["low", "mid", "high"])
out["within_pub_tertile"] = {b: {"n_ess": int(g.ess.sum()), "n_non": int((g.ess == 0).sum()),
                                  "auc_deg700": auc(g.deg700.to_numpy(), g.ess.to_numpy()) if g.ess.nunique() == 2 else None,
                                  "auc_deg700_notm": auc(g.deg700_notm.to_numpy(), g.ess.to_numpy()) if g.ess.nunique() == 2 else None}
                             for b, g in e.groupby("pub_bin", observed=True)}

# fairer test: NEGv1 genes are mostly not expressed in cell lines, which alone gives them few coexpression/experimental
# edges. Restrict to genes broadly expressed in HPA cell lines (nTPM >= 1 in >= 90% of lines) and compare CEGv2 with
# every other broadly expressed gene (that rest still contains some essential genes, so this is conservative).
import zipfile
z = zipfile.ZipFile(f"{D}/rna_celline.tsv.zip")
frac = []
for ch in pd.read_csv(z.open("rna_celline.tsv"), sep="\t", usecols=["Gene name", "Cell line", "nTPM"], chunksize=2_000_000):
    frac.append(ch.assign(on=ch.nTPM >= 1).groupby("Gene name").on.agg(["sum", "count"]))
frac = pd.concat(frac).groupby(level=0).sum()
broad = set(frac.index[(frac["sum"] / frac["count"]) >= 0.9])
out["hpa_broadly_expressed_genes"] = len(broad)
b = m[m.symbol.isin(broad)].copy()
b["ess"] = b.GeneID.isin(set(ceg.ENTREZ_ID)).astype(int)
yb = b.ess.to_numpy()
out["broad_n"] = {"essential": int(yb.sum()), "other": int((yb == 0).sum())}
out["broad_auc_essential"] = {c: auc(b[c].to_numpy(), yb) for c in ("pubs", "deg400", "deg700", "deg400_notm", "deg700_notm")}
out["broad_median"] = {c: {"essential": float(b[b.ess == 1][c].median()), "other": float(b[b.ess == 0][c].median())}
                       for c in ("pubs", "deg700", "deg700_notm")}
lrb = {}
for c in ("deg700", "deg700_notm"):
    X = pd.DataFrame({"degree": np.log1p(b[c]), "pubs": np.log1p(b.pubs)}); X = (X - X.mean()) / X.std()
    fit = sm.Logit(yb, sm.add_constant(X)).fit(disp=0)
    lrb[c] = {"coef_degree": round(float(fit.params["degree"]), 3), "p_degree": float(fit.pvalues["degree"]),
              "coef_pubs": round(float(fit.params["pubs"]), 3), "p_pubs": float(fit.pvalues["pubs"])}
out["broad_logit"] = lrb
b["pub_bin"] = pd.qcut(b.pubs.rank(method="first"), 3, labels=["low", "mid", "high"])
out["broad_within_pub_tertile"] = {k: {"n_ess": int(g.ess.sum()), "n_other": int((g.ess == 0).sum()),
                                       "auc_deg700": auc(g.deg700.to_numpy(), g.ess.to_numpy()),
                                       "auc_deg700_notm": auc(g.deg700_notm.to_numpy(), g.ess.to_numpy())}
                                   for k, g in b.groupby("pub_bin", observed=True)}
# the single-experiment hubs: proteins whose >= 0.4 edges are almost all experimental-only at one score value
e4 = L[L.combined_score >= 400]
expo = e4[(e4[[c for c in CH if c not in ("experimental",)]].sum(axis=1) == 0)]
out["experimental_only_edges_ge400"] = {"n": int(len(expo)), "pct_of_ge400": round(100 * len(expo) / len(e4), 1),
                                        "top_scores": expo.experimental.value_counts().head(5).to_dict()}
def one_value_hub(g):
    sid = deg.loc[deg.symbol == g, "sid"].iloc[0]
    x = e4[(e4.protein1 == sid) | (e4.protein2 == sid)]
    return {"edges_ge400": int(len(x)), "edges_ge700": int((x.combined_score >= 700).sum()),
            "experimental_is_457_pct": round(100 * float((x.experimental == 457).mean()), 1),
            "pubs": float(m.loc[m.symbol == g, "pubs"].iloc[0])}
out["single_experiment_hubs"] = {g: one_value_hub(g) for g in ("TRIM67", "ZRANB1", "PLEKHA4", "TP53")}

# --- review additions ---
# how common is the experimental score 0.457 (is it one dataset's fingerprint?)
x457 = L[L.experimental == 457]
cnt457 = pd.concat([x457.protein1, x457.protein2]).value_counts()
out["exp457"] = {"pairs": int(len(x457)), "proteins": int(len(cnt457)), "proteins_ge100": int((cnt457 >= 100).sum()),
                 "proteins_ge500": int((cnt457 >= 500).sum()),
                 "top10_share_pct": round(100 * float(cnt457.head(10).sum() / (2 * len(x457))), 1),
                 "top": [(deg.set_index("sid").symbol.get(k), int(v)) for k, v in cnt457.head(12).items()]}
# links that need text mining to reach 0.4 but have some other evidence
e4b = L[L.combined_score >= 400]
nzb = e4b[CH] > 0
other = nzb[[c for c in CH if c != "textmining"]].any(axis=1)
out["need_tm_with_other_evidence_pct_ge400"] = round(100 * float(((e4b.no_tm < 400) & other).mean()), 1)
# paper-count AUC within the same publication thirds (broadly expressed genes)
out["broad_within_pub_tertile_auc_pubs"] = {k: auc(g.pubs.to_numpy(), g.ess.to_numpy()) for k, g in b.groupby("pub_bin", observed=True)}
# expression confound: mean nTPM across HPA cell lines
mean_ntpm = []
for ch in pd.read_csv(z.open("rna_celline.tsv"), sep="\t", usecols=["Gene name", "nTPM"], chunksize=2_000_000):
    mean_ntpm.append(ch.groupby("Gene name").nTPM.agg(["sum", "count"]))
mean_ntpm = pd.concat(mean_ntpm).groupby(level=0).sum(); mean_ntpm = mean_ntpm["sum"] / mean_ntpm["count"]
b["expr"] = b.symbol.map(mean_ntpm)
bb = b[b.expr.notna()].copy()
out["broad_auc_expression"] = auc(bb.expr.to_numpy(), bb.ess.to_numpy())
bb["expr_q"] = pd.qcut(bb.expr.rank(method="first"), 5, labels=False)
out["broad_within_expr_quintile_auc_deg700"] = {int(k): auc(g.deg700.to_numpy(), g.ess.to_numpy()) for k, g in bb.groupby("expr_q")}
X = pd.DataFrame({"degree": np.log1p(bb.deg700), "pubs": np.log1p(bb.pubs), "expr": np.log1p(bb.expr)}); X = (X - X.mean()) / X.std()
fit = sm.Logit(bb.ess.to_numpy(), sm.add_constant(X)).fit(disp=0)
out["broad_logit_with_expr"] = {k: round(float(fit.params[k]), 3) for k in ("degree", "pubs", "expr")}
# detected proteome: share broadly expressed, vs all STRING proteins
pg = pd.read_csv("../diann-hye/report.pg_matrix.tsv", sep="\t", usecols=["Protein.Names", "Genes"])
pg = pg[pg["Protein.Names"].str.contains("_HUMAN", na=False) & pg.Genes.notna()]
dgenes = set(pg.Genes.str.split(";").str[0])
dd = deg[deg.symbol.isin(dgenes)]
out["detected_vs_all"] = {"detected_broad_pct": round(100 * float(dd.symbol.isin(broad).mean()), 1),
                          "all_broad_pct": round(100 * float(deg.symbol.isin(broad).mean()), 1),
                          "detected_median_mean_ntpm": float(dd.symbol.map(mean_ntpm).median()),
                          "all_median_mean_ntpm": float(deg.symbol.map(mean_ntpm).median())}
# HPA file date inside the zip
out["hpa_zip_member_date"] = "%04d-%02d-%02d" % z.getinfo("rna_celline.tsv").date_time[:3]
out["hpa_n_cell_lines"] = int(pd.read_csv(z.open("rna_celline.tsv"), sep="\t", usecols=["Cell line"]).iloc[:, 0].nunique())

deg.to_csv("degree_table.tsv.gz", sep="\t", index=False)
json.dump(out, open("hubs_results.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
