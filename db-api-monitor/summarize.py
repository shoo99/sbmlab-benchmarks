"""Summarises checks.csv (written by monitor.py, one row per API per weekly run) into scores.csv; issue #9.

Per database and runner:
  uptime_share       share of runs with HTTP 200
  example_ok_share   share of runs whose response still held the documented example (UniProt returns P00441, ...)
  median_seconds     median latency of the successful runs
  p95_seconds        95th percentile (nearest rank) of the successful runs
  release_changes    how often the reported release changed between runs (the dates are printed)
Also, as the local reference, the one-evening measurement of bioinfo-db-check (2026-10-13, about 23:00 KST, a home
network in Seoul, 2 rounds x 5 calls per API, probe_run1.csv and probe.csv): uptime, median and p95 latency.
GitHub's runners sit in a few cloud regions (mostly the US), so their latency is not what a user in Korea sees.
Usage: python3 summarize.py [checks.csv] [out_dir]   (defaults: this folder; standard library only)
Writes scores.csv into out_dir and prints the tables.
"""
import csv
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from monitor import ENDPOINTS  # noqa: E402  (database names and their order)

checks = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "checks.csv"
out = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE
REFERENCE = [HERE.parent / "bioinfo-db-check" / "probe_run1.csv", HERE.parent / "bioinfo-db-check" / "probe.csv"]
FOLDER = "db-api-monitor"
COLUMNS = ["folder", "dataset", "tool", "version", "setting", "subset", "metric", "value", "n"]
DBS = [e[0] for e in ENDPOINTS]
SUBSET = "main endpoint (human SOD1 example)"


def p95(xs):
    xs = sorted(xs)
    return xs[max(0, -(-95 * len(xs) // 100) - 1)]          # nearest rank


rows = []
weekly = list(csv.DictReader(open(checks))) if checks.exists() else []
runners = sorted({r["runner"] for r in weekly})
for runner in runners:
    w = [r for r in weekly if r["runner"] == runner]
    dates = sorted({r["date_utc"] for r in w})
    dataset = f"weekly checks, {runner}, {dates[0][:10]} to {dates[-1][:10]} ({len(dates)} run{'' if len(dates) == 1 else 's'})"
    print(f"== {dataset}")
    print(f"   {'database':<13} {'runs':>4} {'uptime':>7} {'example':>8} {'median s':>9} {'p95 s':>7}  release (changes)")
    changes_all = []
    for db in DBS:
        r = sorted((x for x in w if x["db"] == db), key=lambda x: x["date_utc"])
        if not r:
            continue
        ok = [float(x["seconds"]) for x in r if x["status"] == "200"]
        rel = [(x["date_utc"][:10], x["release"]) for x in r if x["release"]]
        changes = [(d, a[1], b) for (d, b), a in zip(rel[1:], rel[:-1]) if b != a[1]]
        changes_all += [(db,) + c for c in changes]
        latest = rel[-1][1] if rel else "not reported"
        m = {"uptime_share": round(sum(x["status"] == "200" for x in r) / len(r), 4),
             "example_ok_share": round(sum(x["example_ok"] == "1" for x in r) / len(r), 4)}
        if ok:
            m["median_seconds"] = round(statistics.median(ok), 3)
            m["p95_seconds"] = round(p95(ok), 3)
        if rel:
            m["release_changes"] = len(changes)
        for metric, v in m.items():
            rows.append([FOLDER, dataset, db, f"release {latest}" if rel else latest, "one call per run, weekly",
                         SUBSET, metric, v, len(r)])
        med = f"{m['median_seconds']:.2f}" if ok else "-"
        q = f"{m['p95_seconds']:.2f}" if ok else "-"
        print(f"   {db:<13} {len(r):>4} {m['uptime_share']:>7.3f} {m['example_ok_share']:>8.3f} {med:>9} {q:>7}  "
              f"{latest}{f' ({len(changes)})' if rel else ''}")
    broken = [x for x in w if x["status"] == "200" and x["example_ok"] != "1"]
    failed = [x for x in w if x["status"] != "200"]
    print("   release changes:", "; ".join(f"{db} {a} -> {b} ({d})" for db, d, a, b in changes_all) or "none")
    print("   documented example broken (HTTP 200 but changed content):",
          "; ".join(f"{x['db']} {x['date_utc'][:10]}" for x in broken) or "none")
    print("   failed calls:", "; ".join(f"{x['db']} {x['date_utc'][:10]} {x['status']}" for x in failed) or "none")
if not weekly:
    print("== no weekly checks yet (checks.csv has no rows); only the local reference below")

ref = [r for f in REFERENCE if f.exists() for r in csv.DictReader(open(f))]
dataset = "bioinfo-db-check reference: 2026-10-13 ~23:00 KST, home network in Seoul, 2 rounds x 5 calls"
print(f"\n== {dataset}")
print(f"   {'database':<13} {'calls':>5} {'uptime':>7} {'median s':>9} {'p95 s':>7}")
for db in DBS:
    r = [x for x in ref if x["db"] == db]
    if not r:
        continue
    ok = [float(x["seconds"]) for x in r if x["status"] == "200"]
    m = {"uptime_share": round(len(ok) / len(r), 4)}
    if ok:
        m["median_seconds"] = round(statistics.median(ok), 3)
        m["p95_seconds"] = round(p95(ok), 3)
    for metric, v in m.items():
        rows.append([FOLDER, dataset, db, "as of 2026-10-13", "5 calls 1 s apart, twice", SUBSET, metric, v, len(r)])
    print(f"   {db:<13} {len(r):>5} {m['uptime_share']:>7.3f} {m.get('median_seconds', float('nan')):>9.2f} "
          f"{m.get('p95_seconds', float('nan')):>7.2f}")

out.mkdir(parents=True, exist_ok=True)
with open(out / "scores.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(COLUMNS)
    w.writerows(rows)
print(f"\nwrote scores.csv ({len(rows)} rows)")
