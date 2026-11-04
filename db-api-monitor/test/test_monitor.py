"""Offline test of monitor.py and summarize.py against a local HTTP server (run by CI; standard library only).

The 13 endpoints are pointed at a server on 127.0.0.1 that answers like the real APIs, with each database's own
example and release parsers. Run 1: PRIDE answers HTTP 500. Run 2: UniProt reports a new release and a changed
example. Checked: statuses, latency recorded, every release parsed (header, body or version request), example_ok,
an unreachable host recorded as error:connection, no response body stored in checks.csv, and summarize.py's
uptime, example share and release changes.
Usage: python3 test/test_monitor.py
"""
import csv
import http.server
import json
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import monitor  # noqa: E402

MARK = "BODYMARKER-must-not-be-stored"
state = {"run": 1}
hits = {}


def responses():
    r2 = state["run"] == 2
    return {
        "/uniprot": (200, {"primaryAccession": "P99999" if r2 else "P00441", "x": MARK},
                     {"X-UniProt-Release": "2026_05" if r2 else "2026_04"}),
        "/ensembl": (200, {"id": "ENSG00000142168", "x": MARK}, {}),
        "/ensembl-info": (200, {"releases": [115, 116]}, {}),
        "/gene": (200, {"result": {"uids": ["6647"]}, "x": MARK}, {}),
        "/clinvar": (200, {"esearchresult": {"count": "412"}}, {}),
        "/pdb": (200, {"total_count": 5}, {}),
        "/alphafold": (200, [{"entryId": "AF-P00441-F1", "latestVersion": 6}], {}),
        "/interpro": (200, {"count": 7}, {}),
        "/interpro-root": (200, {"databases": {"interpro": {"version": "105.0"}}}, {}),
        "/quickgo": (200, {"numberOfHits": 300}, {}),
        "/string": (200, [{"preferredName_A": "SOD1", "x": MARK}], {}),
        # run 1: the STRING version request fails once and succeeds on the retry
        "/string-version": ((503, {"e": 1}, {}) if not r2 and hits.get("/string-version", 0) == 1
                            else (200, [{"string_version": "12.5"}], {})),
        "/reactome": (200, [{"stId": "R-HSA-1"}], {}),
        # run 2: Reactome's version request answers 200 with nothing readable (a parser problem)
        "/reactome-version": (200, "" if r2 else "97", {}),
        "/kegg": (200, "hsa:6647\tpath:hsa04213\n", {}),
        # run 1: KEGG's info request fails on both tries (a path or server problem); run 2: the real format
        "/kegg-info": ((403, {"e": 1}, {}) if not r2 else
                       (200, "pathway          KEGG Pathway Database\npath             Release 118.0+/04-08, Apr 26\n"
                             "                 Kanehisa Laboratories\n                 585 entries\n", {})),
        "/geo": (200, {"esearchresult": {"count": "30"}}, {}),
        "/pride": (200, [{"accession": "PXD000001"}], {}) if r2 else (500, {"error": MARK}, {}),
    }


class Handler(http.server.BaseHTTPRequestHandler):
    def _answer(self):
        hits[self.path] = hits.get(self.path, 0) + 1
        status, body, headers = responses()[self.path]
        data = body.encode() if isinstance(body, str) else json.dumps(body).encode()
        self.send_response(status)
        for k, v in headers.items():
            self.send_header(k, v)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    do_GET = do_POST = _answer

    def log_message(self, *a):
        pass


server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
threading.Thread(target=server.serve_forever, daemon=True).start()
base = f"http://127.0.0.1:{server.server_port}"
paths = {"UniProt": ("/uniprot", None), "Ensembl": ("/ensembl", "/ensembl-info"), "NCBI Gene": ("/gene", None),
         "ClinVar": ("/clinvar", None), "RCSB PDB": ("/pdb", None), "AlphaFold DB": ("/alphafold", None),
         "InterPro": ("/interpro", "/interpro-root"), "QuickGO": ("/quickgo", None),
         "STRING": ("/string", "/string-version"), "Reactome": ("/reactome", "/reactome-version"),
         "KEGG": ("/kegg", "/kegg-info"), "GEO": ("/geo", None), "PRIDE": ("/pride", None)}
eps = []
for name, method, url, body, example, release in monitor.ENDPOINTS:
    main_path, version_path = paths[name]
    if release and release[0] not in ("header", "body"):
        release = (base + version_path, release[1])
    eps.append((name, method, base + main_path, body, example, release))
monitor.ENDPOINTS = eps
monitor.time.sleep = lambda s: None                      # no 1 s pauses in the test
checks = []
with tempfile.TemporaryDirectory() as tmp:
    t = Path(tmp)
    out = t / "checks.csv"
    for run in (1, 2):
        state["run"] = run
        sys.argv = ["monitor.py", str(out)]
        monitor.main()
        time.sleep(0.01)
    rows = list(csv.DictReader(open(out)))
    r1 = {r["db"]: r for r in rows[:13]}
    r2 = {r["db"]: r for r in rows[13:]}
    checks.append(("26 rows, header written once", len(rows) == 26 and open(out).read().count("date_utc") == 1))
    want = {"UniProt": "2026_04", "Ensembl": "116", "AlphaFold DB": "6", "InterPro": "105.0", "STRING": "12.5",
            "Reactome": "97", "KEGG": "", "NCBI Gene": "", "PRIDE": ""}
    for db, rel in want.items():
        checks.append((f"run 1 {db}: release '{rel}'", r1[db]["release"] == rel))
    checks.append(("run 1 release_status: header / body / version request status",
                   (r1["UniProt"]["release_status"], r1["AlphaFold DB"]["release_status"], r1["Ensembl"]["release_status"],
                    r1["NCBI Gene"]["release_status"]) == ("header", "body", "200", "")))
    checks.append(("run 1 STRING: version request failed once, retry succeeded (503>200)",
                   r1["STRING"]["release_status"] == "503>200"))
    checks.append(("run 1 KEGG: version request failed twice (403>403), release empty",
                   r1["KEGG"]["release_status"] == "403>403" and r1["KEGG"]["release"] == ""))
    checks.append(("run 2 KEGG: real info format parsed as 118.0+", r2["KEGG"]["release"] == "118.0+"))
    checks.append(("run 2 Reactome: answered 200 but nothing readable", r2["Reactome"]["release_status"] == "200"
                   and r2["Reactome"]["release"] == ""))
    checks.append(("run 1: twelve APIs HTTP 200 with the example intact",
                   sum(r["status"] == "200" and r["example_ok"] == "1" for r in r1.values()) == 12))
    checks.append(("run 1 PRIDE: status 500, example_ok 0", r1["PRIDE"]["status"] == "500" and r1["PRIDE"]["example_ok"] == "0"))
    checks.append(("run 2 UniProt: new release 2026_05, example broken (status 200, example_ok 0)",
                   r2["UniProt"]["release"] == "2026_05" and r2["UniProt"]["status"] == "200" and r2["UniProt"]["example_ok"] == "0"))
    checks.append(("latency recorded for every call", all(float(r["seconds"]) >= 0 for r in rows)))
    checks.append(("no response body stored", MARK not in open(out).read()))
    port = server.server_port
    server.shutdown()
    server.server_close()
    status, sec, _, _ = monitor.call("GET", f"http://127.0.0.1:{port}/uniprot", None, timeout=2)
    checks.append((f"unreachable host recorded as error:connection ({status})", status == "error:connection"))

    summ = subprocess.run([sys.executable, str(HERE / "summarize.py"), str(out), str(t)], check=True, capture_output=True,
                          text=True).stdout
    checks.append(("summary lists the failed KEGG version request", "version request failed: KEGG" in summ))
    checks.append(("summary lists the STRING retry", "needed a retry: STRING" in summ))
    checks.append(("summary lists Reactome as answered but not read", "no release read: Reactome" in summ))
    sc = list(csv.DictReader(open(t / "scores.csv")))
    v = lambda db, metric: next(float(r["value"]) for r in sc if r["tool"] == db and r["metric"] == metric
                                and r["dataset"].startswith("weekly"))
    checks.append(("summary: PRIDE uptime 0.5", v("PRIDE", "uptime_share") == 0.5))
    checks.append(("summary: UniProt example_ok share 0.5, one release change",
                   v("UniProt", "example_ok_share") == 0.5 and v("UniProt", "release_changes") == 1))
    checks.append(("summary: KEGG release seen once, no change counted", v("KEGG", "release_changes") == 0))
    checks.append(("summary: Ensembl uptime 1, no release change",
                   v("Ensembl", "uptime_share") == 1 and v("Ensembl", "release_changes") == 0))
    checks.append(("summary: reference rows from bioinfo-db-check present",
                   any(r["dataset"].startswith("bioinfo-db-check reference") for r in sc)))

for name, ok in checks:
    print(("ok   " if ok else "FAIL ") + name)
if not all(ok for _, ok in checks):
    sys.exit(1)
print(f"all {len(checks)} checks passed")
