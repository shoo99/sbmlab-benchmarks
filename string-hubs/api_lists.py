"""What a user sees when they paste a protein list into STRING (12.5), using proteins actually detected by mass spec.

Detected set: human protein groups in ../diann-hye/report.pg_matrix.tsv (first gene of each group), mapped to STRING 12.5 by
preferred name (local protein.info). Lists:
  random: 20 draws of 200 detected proteins (seeds 1-20)
  planted: 20 draws of 15 detected Reactome 'Cholesterol biosynthesis' (R-HSA-191273, Reactome v97) proteins + 185 random
For each list: STRING API ppi_enrichment at required_score 400 and 700, without and with the detected set as background;
and, from the local network file, the 10 proteins with most partners inside the list (the 'hubs' a user would look at).
Outputs: api_results.json
"""
import json, time, random, urllib.request, urllib.parse, sys
import numpy as np, pandas as pd

D = "/home/sysoft/data/string125"
API = "https://version-12-5.string-db.org/api/json/ppi_enrichment"
CHOL = ['ACAT2', 'ARV1', 'CYP51A1', 'DHCR24', 'DHCR7', 'EBP', 'FDFT1', 'FDPS', 'GGPS1', 'HMGCR', 'HMGCS1', 'HSD17B7', 'IDI1',
        'IDI2', 'LBR', 'LSS', 'MSMO1', 'MVD', 'MVK', 'NSDHL', 'PLPP6', 'PMVK', 'SC5D', 'SQLE', 'SREBF1', 'SREBF2', 'TM7SF2']

pg = pd.read_csv("../diann-hye/report.pg_matrix.tsv", sep="\t", usecols=["Protein.Names", "Genes"])
pg = pg[pg["Protein.Names"].str.contains("_HUMAN", na=False) & pg.Genes.notna()]
genes = pg.Genes.str.split(";").str[0].unique()
info = pd.read_csv(f"{D}/9606.protein.info.v12.5.txt.gz", sep="\t").rename(columns={"#string_protein_id": "sid"})
name2sid = dict(zip(info.preferred_name, info.sid))
det = sorted({name2sid[g] for g in genes if g in name2sid})
chol = sorted({name2sid[g] for g in CHOL if g in name2sid and name2sid[g] in set(det)})
other = sorted(set(det) - set(chol))
sid2name = dict(zip(info.sid, info.preferred_name))

deg = pd.read_csv("degree_table.tsv.gz", sep="\t")
pubs = dict(zip(deg.sid, deg.pubs))
L = pd.read_csv(f"{D}/9606.protein.links.detailed.v12.5.txt.gz", sep=" ", usecols=["protein1", "protein2", "combined_score"])
L = L[L.protein1 < L.protein2]
detset = set(det)
L = L[L.protein1.isin(detset) & L.protein2.isin(detset)]

def call(ids, score, background):
    data = {"identifiers": "%0d".join(ids), "species": 9606, "required_score": score, "caller_identity": "sbmlab_benchmark"}
    if background:
        data["background_string_identifiers"] = "%0d".join(det)
    body = urllib.parse.urlencode(data).encode()
    for attempt in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(API, data=body), timeout=120) as r:
                return json.load(r)[0]
        except Exception as e:
            time.sleep(5 * (attempt + 1)); err = e
    raise err

def local_hubs(ids, score):
    s = set(ids)
    e = L[(L.combined_score >= score) & L.protein1.isin(s) & L.protein2.isin(s)]
    d = pd.concat([e.protein1, e.protein2]).value_counts()
    return [(sid2name[k], int(v), pubs.get(k)) for k, v in d.head(10).items()]

res = {"n_detected_genes": int(len(genes)), "n_detected_string": len(det), "n_chol_detected": len(chol),
       "chol_detected": [sid2name[c] for c in chol], "random": [], "planted": []}
for kind in ("random", "planted"):
    for seed in range(1, 21):
        rng = random.Random(seed)
        ids = rng.sample(det, 200) if kind == "random" else rng.sample(chol, 15) + rng.sample(other, 185)
        row = {"seed": seed}
        for score in (400, 700):
            for bg in (False, True):
                r = call(ids, score, bg)
                row[f"s{score}_{'bg' if bg else 'genome'}"] = {"edges": r["number_of_edges"], "expected": r["expected_number_of_edges"],
                                                              "p": r["p_value"]}
                time.sleep(1)
            h = local_hubs(ids, score)
            row[f"hubs{score}"] = h
            if kind == "planted":
                cs = {sid2name[c] for c in chol}
                row[f"chol_in_top10_{score}"] = sum(1 for n, _, _ in h if n in cs)
        row["median_pubs_list"] = float(np.median([pubs.get(i, 0) or 0 for i in ids]))
        res[kind].append(row)
        print(kind, seed, {k: v for k, v in row.items() if k.startswith("s")}, row.get("chol_in_top10_400"), row.get("chol_in_top10_700"), flush=True)

def summ(rows, key):
    p = np.array([r[key]["p"] for r in rows]); ratio = np.array([r[key]["edges"] / max(r[key]["expected"], 1e-9) for r in rows])
    return {"n_p_lt_0.05": int((p < 0.05).sum()), "median_p": float(np.median(p)), "median_obs_over_exp": round(float(np.median(ratio)), 2)}
res["summary"] = {k: {key: summ(res[k], key) for key in ("s400_genome", "s400_bg", "s700_genome", "s700_bg")} for k in ("random", "planted")}
res["summary"]["planted_chol_in_top10"] = {s: [r[f"chol_in_top10_{s}"] for r in res["planted"]] for s in (400, 700)}
hp = [h[2] for r in res["random"] for h in r["hubs400"][:3] if h[2] is not None]
res["summary"]["random_top3_hub_median_pubs_400"] = float(np.median(hp))
res["summary"]["random_list_median_pubs"] = float(np.median([r["median_pubs_list"] for r in res["random"]]))
json.dump(res, open("api_results.json", "w"), indent=1)
print(json.dumps(res["summary"], indent=1))
