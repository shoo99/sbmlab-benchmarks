"""Repair the publication-count fields of the bundled api_results.json without calling the STRING API (2026-11-02).

api_lists.py used `pubs.get(i, 0) or 0`, which lets NaN through (NaN is truthy) for proteins without an NCBI GeneID,
so 14 of the 40 lists got median_pubs_list = NaN, summary.random_list_median_pubs was NaN, one in-list hub (MT-CYB,
planted seed 7) had pubs NaN, and the file held bare NaN, which strict JSON parsers reject. This script rebuilds the
same 40 protein lists (same detected set and seeds as api_lists.py, with STRING ids and preferred names taken from
degree_table.tsv.gz instead of the STRING info file), checks that the 26 medians that were not NaN come out
unchanged, then skips proteins without a GeneID (as hubs.py does), writes NaN hub counts as null and rewrites
api_results.json as strict JSON. API results (edges, expected, p) are untouched.

Usage: cd string-hubs && python3 repair_pubs.py
"""
import json, math, random
import numpy as np, pandas as pd

CHOL = ['ACAT2', 'ARV1', 'CYP51A1', 'DHCR24', 'DHCR7', 'EBP', 'FDFT1', 'FDPS', 'GGPS1', 'HMGCR', 'HMGCS1', 'HSD17B7', 'IDI1',
        'IDI2', 'LBR', 'LSS', 'MSMO1', 'MVD', 'MVK', 'NSDHL', 'PLPP6', 'PMVK', 'SC5D', 'SQLE', 'SREBF1', 'SREBF2', 'TM7SF2']
pg = pd.read_csv("../diann-hye/report.pg_matrix.tsv", sep="\t", usecols=["Protein.Names", "Genes"])
pg = pg[pg["Protein.Names"].str.contains("_HUMAN", na=False) & pg.Genes.notna()]
genes = pg.Genes.str.split(";").str[0].unique()
deg = pd.read_csv("degree_table.tsv.gz", sep="\t")
name2sid = dict(zip(deg.symbol, deg.sid))
det = sorted({name2sid[g] for g in genes if g in name2sid})
chol = sorted({name2sid[g] for g in CHOL if g in name2sid and name2sid[g] in set(det)})
other = sorted(set(det) - set(chol))
pubs = dict(zip(deg.sid, deg.pubs))


def parse_nan(x):
    return float("nan")


res = json.load(open("api_results.json"), parse_constant=parse_nan)
assert (res["n_detected_string"], res["n_chol_detected"]) == (len(det), len(chol)), "detected set differs from api_lists.py"
fixed = 0
for kind in ("random", "planted"):
    for row in res[kind]:
        rng = random.Random(row["seed"])
        ids = rng.sample(det, 200) if kind == "random" else rng.sample(chol, 15) + rng.sample(other, 185)
        known = [pubs[i] for i in ids if pubs[i] == pubs[i]]
        new = float(np.median(known))
        old = row["median_pubs_list"]
        if old == old:
            assert new == old, (kind, row["seed"], old, new)   # lists rebuilt exactly
        else:
            row["median_pubs_list"] = new
            fixed += 1
            print(f"{kind} seed {row['seed']}: median_pubs_list NaN -> {new} ({200 - len(known)} of 200 without GeneID)")
        for s in (400, 700):
            for h in row[f"hubs{s}"]:
                if isinstance(h[2], float) and math.isnan(h[2]):
                    print(f"{kind} seed {row['seed']}: hub {h[0]} pubs NaN -> null")
                    h[2] = None
res["summary"]["random_list_median_pubs"] = float(np.median([r["median_pubs_list"] for r in res["random"]]))
print("lists repaired:", fixed, "| summary.random_list_median_pubs:", res["summary"]["random_list_median_pubs"])
json.dump(res, open("api_results.json", "w"), indent=1, allow_nan=False)
