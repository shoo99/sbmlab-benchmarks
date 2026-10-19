"""Extra summary of api_results.json (no API calls): publication counts of in-list hubs vs the detected set."""
import json, numpy as np, pandas as pd
r = json.load(open("api_results.json"))
deg = pd.read_csv("degree_table.tsv.gz", sep="\t")
pg = pd.read_csv("../diann-hye/report.pg_matrix.tsv", sep="\t", usecols=["Protein.Names", "Genes"])
pg = pg[pg["Protein.Names"].str.contains("_HUMAN", na=False) & pg.Genes.notna()]
genes = set(pg.Genes.str.split(";").str[0])
det = deg[deg.symbol.isin(genes)]
out = {"detected_median_pubs": float(det.pubs.median()), "all_string_median_pubs": float(deg.pubs.median())}
for s in (400, 700):
    top3 = [h[2] for x in r["random"] for h in x[f"hubs{s}"][:3] if h[2] == h[2] and h[2] is not None]
    out[f"random_top3_hub_median_pubs_{s}"] = float(np.median(top3))
    ch = r["summary"]["planted_chol_in_top10"][str(s)]
    out[f"planted_chol_top10_{s}"] = {"median": float(np.median(ch)), "min": int(min(ch)), "max": int(max(ch))}
json.dump(out, open("api_summary_extra.json", "w"), indent=1); print(out)
