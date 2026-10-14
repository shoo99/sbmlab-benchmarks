"""Follow one protein (human SOD1, UniProt P00441, NCBI Gene 6647) through public bioinformatics database APIs.
Each endpoint is called N times, 1 s apart; we record HTTP status and wall-clock latency. Results -> probe.csv / probe.json."""
import json, time, statistics, urllib.request, urllib.parse, csv, datetime, sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
UA = {"User-Agent": "sbmlab-db-check/1.0 (blog benchmark; contact via sbmlab.com)"}
PDB_Q = {"query": {"type": "terminal", "service": "text", "parameters": {
    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
    "operator": "exact_match", "value": "P00441"}}, "return_type": "entry", "request_options": {"return_all_hits": True}}
E = [  # name, method, url, body
    ("UniProt", "GET", "https://rest.uniprot.org/uniprotkb/P00441.json", None),
    ("Ensembl", "GET", "https://rest.ensembl.org/lookup/symbol/homo_sapiens/SOD1?content-type=application/json", None),
    ("NCBI Gene", "GET", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=6647&retmode=json", None),
    ("ClinVar", "GET", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=SOD1%5Bgene%5D&retmode=json", None),
    ("RCSB PDB", "POST", "https://search.rcsb.org/rcsbsearch/v2/query", PDB_Q),
    ("AlphaFold DB", "GET", "https://alphafold.ebi.ac.uk/api/prediction/P00441", None),
    ("InterPro", "GET", "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P00441", None),
    ("QuickGO", "GET", "https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId=P00441&limit=1", None),
    ("STRING", "GET", "https://string-db.org/api/json/interaction_partners?identifiers=SOD1&species=9606&limit=10", None),
    ("Reactome", "GET", "https://reactome.org/ContentService/data/mapping/UniProt/P00441/pathways?species=9606", None),
    ("KEGG", "GET", "https://rest.kegg.jp/link/pathway/hsa:6647", None),
    ("GEO", "GET", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=SOD1%20AND%20%22Homo%20sapiens%22%5Borgn%5D&retmode=json", None),
    ("PRIDE", "GET", "https://www.ebi.ac.uk/pride/ws/archive/v3/search/projects?keyword=SOD1&pageSize=100", None),
]

def call(method, url, body):
    data = json.dumps(body).encode() if body else None
    req = urllib.request.Request(url, data=data, method=method, headers={**UA, **({"Content-Type": "application/json"} if body else {})})
    t = time.perf_counter()
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            payload = r.read(); status = r.status; hdr = dict(r.headers)
    except urllib.error.HTTPError as e:
        payload, status, hdr = e.read(), e.code, dict(e.headers)
    except Exception as e:
        return None, round(time.perf_counter() - t, 3), None, str(e)[:80]
    return status, round(time.perf_counter() - t, 3), payload, hdr

rows, last = [], {}
for name, m, url, body in E:
    for i in range(N):
        status, sec, payload, hdr = call(m, url, body)
        rows.append({"db": name, "try": i + 1, "status": status, "seconds": sec})
        if status == 200: last[name] = (payload, hdr)
        time.sleep(1)
    ok = [r for r in rows if r["db"] == name and r["status"] == 200]
    print(f"{name:13s} ok {len(ok)}/{N}  median {statistics.median([r['seconds'] for r in ok]) if ok else float('nan'):.2f} s", flush=True)

stamp = datetime.datetime.now().astimezone().isoformat(timespec="minutes")
with open("probe.csv", "w", newline="") as f:
    w = csv.DictWriter(f, ["db", "try", "status", "seconds"]); w.writeheader(); w.writerows(rows)
json.dump({k: {"body": v[0].decode("utf8", "replace"), "headers": {h: v[1][h] for h in v[1] if h.lower().startswith(("x-", "content-type", "last-modified"))}} for k, v in last.items()},
          open("probe.json", "w"), indent=1)
print("run at", stamp)
