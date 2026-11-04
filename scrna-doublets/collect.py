"""Collects the per-droplet truth and tool scores of every prepared data set into two small tables for bundling.

  <out_dir>/cell_scores.tsv.gz   one row per droplet: dataset, barcode, truth, then <tool> score and <tool> call
                                 columns (call: 1 / 0, empty when the tool makes no call); gzip without a time stamp
  <out_dir>/datasets.tsv         info.tsv of each data set, plus each tool's version and setting
score.py scores these two files, so the scoring re-runs from the repository alone.
Usage: python3 collect.py <work_dir> [out_dir]    (out_dir defaults to this folder; standard library only)
"""
import csv
import gzip
import io
import sys
from pathlib import Path

work = Path(sys.argv[1])
out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).resolve().parent
sets = sorted(p.parent for p in work.glob("*/truth.tsv"))
if not sets:
    sys.exit(f"no prepared data sets (*/truth.tsv) under {work}")

info_rows, tools_all, cells = [], [], []
for d in sets:
    info = next(csv.DictReader(open(d / "info.tsv"), delimiter="\t"))
    truth = {r["barcode"]: r["truth"] for r in csv.DictReader(open(d / "truth.tsv"), delimiter="\t")}
    scores, versions = {}, {}
    for f in sorted(d.glob("scores_*.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            scores.setdefault(r["tool"], {})[r["barcode"]] = (r["score"], {"TRUE": "1", "FALSE": "0"}.get(r["call"], ""))
    for f in sorted(d.glob("versions_*.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            versions[r["tool"]] = (r["version"], r.get("setting", ""))
    for t in scores:
        if t not in tools_all:
            tools_all.append(t)
        missing = set(truth) - set(scores[t])
        if missing:
            sys.exit(f"{d.name}: {t} has no score for {len(missing)} droplets")
    info["tools"] = "; ".join(f"{t} = {versions.get(t, ('?', ''))[0]}" + (f" ({versions[t][1]})" if versions.get(t, ('', ''))[1] else "")
                              for t in scores)
    info_rows.append(info)
    cells.append((d.name, truth, scores))

out.mkdir(parents=True, exist_ok=True)
with open(out / "datasets.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(info_rows[0].keys()), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(info_rows)
buf = io.StringIO()
w = csv.writer(buf, delimiter="\t", lineterminator="\n")
w.writerow(["dataset", "barcode", "truth"] + [f"{t} {k}" for t in tools_all for k in ("score", "call")])
n = 0
for name, truth, scores in cells:
    for b, tr in truth.items():
        row = [name, b, tr]
        for t in tools_all:
            row += list(scores[t][b]) if t in scores else ["", ""]
        w.writerow(row)
        n += 1
with open(out / "cell_scores.tsv.gz", "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz:
    gz.write(buf.getvalue().encode())
print(f"wrote cell_scores.tsv.gz ({n} droplets, {len(cells)} data sets, {len(tools_all)} tools) and datasets.tsv")
