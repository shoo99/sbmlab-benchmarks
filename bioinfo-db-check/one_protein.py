import json, time, http.client, urllib.request, urllib.error

def get(url, body=None, tries=4):
    """GET (or POST JSON) with retries. Retry only what can recover: 429, 5xx and network errors."""
    for i in range(tries):
        headers = {"User-Agent": "db-walkthrough/1.0"}
        if body:
            headers["Content-Type"] = "application/json"
        req = urllib.request.Request(url, data=json.dumps(body).encode() if body else None, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return r.read().decode()
        except urllib.error.HTTPError as err:
            if err.code != 429 and err.code < 500:
                raise                              # 400/404: a bad ID or query, retrying won't help
            wait = float(err.headers.get("Retry-After") or 2 ** i)
        except (urllib.error.URLError, OSError, http.client.HTTPException):
            wait = 2 ** i                          # timeouts, resets, dropped connections
        if i == tries - 1:
            raise RuntimeError(f"gave up after {tries} tries: {url}")
        time.sleep(wait)                           # 1, 2, 4 s unless the server says otherwise

acc = "P00441"                                      # human SOD1

u = json.loads(get(f"https://rest.uniprot.org/uniprotkb/{acc}.json"))
print("UniProt   ", u["uniProtkbId"], u["proteinDescription"]["recommendedName"]["fullName"]["value"],
      "|", u["sequence"]["length"], "aa")

e = json.loads(get("https://rest.ensembl.org/lookup/symbol/homo_sapiens/SOD1?content-type=application/json"))
print("Ensembl   ", e["id"], f"chr{e['seq_region_name']}:{e['start']}-{e['end']}", e["assembly_name"])

pdb = json.loads(get("https://search.rcsb.org/rcsbsearch/v2/query", {
    "query": {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        "operator": "exact_match", "value": acc}},
    "return_type": "entry", "request_options": {"return_all_hits": True}}))
print("RCSB PDB  ", pdb["total_count"], "entries, e.g.", ", ".join(h["identifier"] for h in pdb["result_set"][:3]))

af = json.loads(get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"))[0]
print("AlphaFold ", af["entryId"], "v%s" % af["latestVersion"], "| mean pLDDT", af["globalMetricValue"])

s = json.loads(get("https://string-db.org/api/json/interaction_partners?identifiers=SOD1&species=9606&limit=5"))
print("STRING    ", ", ".join(f"{x['preferredName_B']} {x['score']}" for x in s))

r = json.loads(get(f"https://reactome.org/ContentService/data/mapping/UniProt/{acc}/pathways?species=9606"))
print("Reactome  ", "; ".join(p["displayName"].strip() for p in r))

c = json.loads(get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=SOD1[gene]&retmode=json"))
print("ClinVar   ", c["esearchresult"]["count"], "records for SOD1")
