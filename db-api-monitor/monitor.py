"""One light check of 13 public bioinformatics database APIs; appends one row per API to checks.csv (issue #9).

Same endpoints as bioinfo-db-check/probe.py (human SOD1: UniProt P00441, NCBI Gene 6647), each called ONCE per run,
1 s apart, plus one small version request for the five databases that report a release that way. Per API we store:
  date_utc    start of the run (UTC, minutes)
  runner      where it ran (MONITOR_RUNNER, e.g. "github-actions ubuntu-24.04"; default "local")
  db          database
  status      HTTP status, or error:<kind> (timeout, dns, connection, ...)
  seconds     wall-clock time of the main request
  release     release / version the API reports (empty when it reports none or the version request failed)
  release_status  where the release came from: "header" or "body" of the main response, or the HTTP status of the
              version request (a failed first try and its retry 3 s later are both shown, e.g. "403>200"); empty
              when the API has no release source or the main request failed
  example_ok  1 if the response still holds the documented example (e.g. UniProt returns P00441), else 0
Response bodies are never stored (some sources, e.g. KEGG, do not allow redistribution); only the fields above.
User-Agent names the project, without a personal contact.
Usage: python3 monitor.py [checks.csv]     (default: checks.csv next to this script; standard library only)
"""
import csv
import datetime
import json
import os
import re
import socket
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

UA = "sbmlab-benchmarks-db-api-monitor/1.0 (+https://github.com/shoo99/sbmlab-benchmarks)"
COLUMNS = ["date_utc", "runner", "db", "status", "seconds", "release", "release_status", "example_ok"]
PDB_Q = {"query": {"type": "terminal", "service": "text", "parameters": {
    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
    "operator": "exact_match", "value": "P00441"}}, "return_type": "entry", "request_options": {"return_all_hits": True}}


def jget(body):
    return json.loads(body.decode("utf8", "replace"))


def kegg_release(body):
    m = re.search(r"Release\s+(\S+?)(?:/|,|\s|$)", body.decode("utf8", "replace"))
    return m.group(1) if m else ""


# name, method, url, body, example check on the main response, release: ("header", name) | ("body", fn) | (url, fn)
ENDPOINTS = [
    ("UniProt", "GET", "https://rest.uniprot.org/uniprotkb/P00441.json", None,
     lambda b: jget(b)["primaryAccession"] == "P00441", ("header", "X-UniProt-Release")),
    ("Ensembl", "GET", "https://rest.ensembl.org/lookup/symbol/homo_sapiens/SOD1?content-type=application/json", None,
     lambda b: jget(b)["id"] == "ENSG00000142168",
     ("https://rest.ensembl.org/info/data?content-type=application/json", lambda b: str(max(jget(b)["releases"])))),
    ("NCBI Gene", "GET", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=6647&retmode=json", None,
     lambda b: "6647" in jget(b)["result"]["uids"], None),
    ("ClinVar", "GET", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=SOD1%5Bgene%5D&retmode=json",
     None, lambda b: int(jget(b)["esearchresult"]["count"]) > 0, None),
    ("RCSB PDB", "POST", "https://search.rcsb.org/rcsbsearch/v2/query", PDB_Q,
     lambda b: jget(b)["total_count"] > 0, None),
    ("AlphaFold DB", "GET", "https://alphafold.ebi.ac.uk/api/prediction/P00441", None,
     lambda b: jget(b)[0]["entryId"].startswith("AF-P00441"), ("body", lambda b: str(jget(b)[0]["latestVersion"]))),
    ("InterPro", "GET", "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P00441", None,
     lambda b: jget(b)["count"] > 0,
     ("https://www.ebi.ac.uk/interpro/api/", lambda b: str(jget(b)["databases"]["interpro"]["version"]))),
    ("QuickGO", "GET", "https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId=P00441&limit=1", None,
     lambda b: jget(b)["numberOfHits"] > 0, None),
    ("STRING", "GET", "https://string-db.org/api/json/interaction_partners?identifiers=SOD1&species=9606&limit=10", None,
     lambda b: len(jget(b)) > 0 and "preferredName_A" in jget(b)[0],
     ("https://string-db.org/api/json/version", lambda b: str(jget(b)[0]["string_version"]))),
    ("Reactome", "GET", "https://reactome.org/ContentService/data/mapping/UniProt/P00441/pathways?species=9606", None,
     lambda b: len(jget(b)) > 0 and "stId" in jget(b)[0],
     ("https://reactome.org/ContentService/data/database/version", lambda b: b.decode().strip())),
    ("KEGG", "GET", "https://rest.kegg.jp/link/pathway/hsa:6647", None,
     lambda b: b"path:hsa" in b, ("https://rest.kegg.jp/info/pathway", kegg_release)),
    ("GEO", "GET", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=SOD1%20AND%20%22Homo%20sapiens%22%5Borgn%5D&retmode=json",
     None, lambda b: int(jget(b)["esearchresult"]["count"]) > 0, None),
    ("PRIDE", "GET", "https://www.ebi.ac.uk/pride/ws/archive/v3/search/projects?keyword=SOD1&pageSize=100", None,
     lambda b: len(jget(b)) > 0, None),
]


def error_kind(e):
    if isinstance(e, (socket.timeout, TimeoutError)) or "timed out" in str(e):
        return "error:timeout"
    reason = getattr(e, "reason", e)
    if isinstance(reason, socket.gaierror):
        return "error:dns"
    if isinstance(reason, (ConnectionError, OSError)):
        return "error:connection"
    return f"error:{type(e).__name__}"


def call(method, url, body, timeout=30):
    """Returns (status, seconds, payload, headers); payload and headers are None after a network error."""
    data = json.dumps(body).encode() if body else None
    headers = {"User-Agent": UA, **({"Content-Type": "application/json"} if body else {})}
    req = urllib.request.Request(url, data=data, method=method, headers=headers)
    t = time.perf_counter()
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            payload, status, hdr = r.read(), str(r.status), r.headers
    except urllib.error.HTTPError as e:
        payload, status, hdr = e.read(), str(e.code), e.headers
    except Exception as e:  # noqa: BLE001 - every failure is a result here
        return error_kind(e), round(time.perf_counter() - t, 3), None, None
    return status, round(time.perf_counter() - t, 3), payload, hdr


def check(endpoint, pause=1.0):
    """Returns (status, seconds, release, release_status, example_ok) for one API."""
    name, method, url, body, example, release = endpoint
    status, sec, payload, hdr = call(method, url, body)
    ok = 0
    if status == "200":
        try:
            ok = int(bool(example(payload)))
        except Exception:  # noqa: BLE001 - a changed response shape is what example_ok reports
            ok = 0
    rel, rel_status = "", ""
    if release and status == "200":
        kind, how = release
        try:
            if kind == "header":
                rel_status = "header"
                rel = hdr.get(how, "") or ""
            elif kind == "body":
                rel_status = "body"
                rel = how(payload)
            else:
                time.sleep(pause)
                s2, _, p2, _ = call("GET", kind, None)
                tries = [s2]
                if s2 != "200":                      # one retry, 3 s later (rate limits, transient errors)
                    time.sleep(3 * pause)
                    s2, _, p2, _ = call("GET", kind, None)
                    tries.append(s2)
                rel_status = ">".join(tries)
                rel = how(p2) if s2 == "200" else ""
        except Exception as e:  # noqa: BLE001
            rel, rel_status = "", rel_status or f"error:{type(e).__name__}"
    return status, sec, rel.strip()[:40], rel_status, ok


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "checks.csv"
    runner = os.environ.get("MONITOR_RUNNER", "local")
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%MZ")
    rows = []
    for ep in ENDPOINTS:
        status, sec, rel, rel_status, ok = check(ep)
        rows.append([stamp, runner, ep[0], status, f"{sec:.3f}", rel, rel_status, ok])
        print(f"{ep[0]:<13} {status:<17} {sec:>7.2f} s  release {rel or '-':<12} ({rel_status or 'no source'})  "
              f"example {'ok' if ok else 'FAILED'}", flush=True)
        time.sleep(1)
    new = not out.exists() or out.stat().st_size == 0
    with open(out, "a", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        if new:
            w.writerow(COLUMNS)
        w.writerows(rows)
    print(f"appended {len(rows)} rows to {out.name} ({stamp}, runner {runner})")


if __name__ == "__main__":
    main()
