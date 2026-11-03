"""Start a new analysis folder from docs/template and register it in catalog.csv (standard library only).

Usage: python3 ci/new_folder.py <folder> --area <area> --post <url> [--title "Question the folder answers"]
  folder  lowercase name like rnaseq-quant-seqc (letters, digits, '-', '.'); must not exist yet
  area    e.g. proteomics, transcriptomics, single-cell, statistics, microbiome, networks-enrichment, databases
  post    the post, or the tracking issue until the post exists (https://...)

Creates <folder>/README.md (template, title filled in) and <folder>/run.sh, then adds a catalog.csv row with
ci = no and TODO in data, truth and tools. ci/lint.py rejects TODO, so CI stays red until the row is filled in.
Afterwards: python3 docs/make_summaries.py to refresh docs/catalog.md.
"""
import argparse
import csv
import re
import shutil
import stat
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
COLUMNS = ["folder", "area", "data", "truth", "tools", "post", "ci"]

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("folder")
ap.add_argument("--area", required=True)
ap.add_argument("--post", required=True)
ap.add_argument("--title", default=None)
a = ap.parse_args()

if not re.fullmatch(r"[a-z0-9][a-z0-9.-]*", a.folder):
    sys.exit(f"folder name '{a.folder}': use lowercase letters, digits, '-' and '.'")
if (ROOT / a.folder).exists():
    sys.exit(f"{a.folder}/ already exists")
if a.folder in {"docs", "ci", "common"}:
    sys.exit(f"{a.folder} is reserved")
if not a.post.startswith("https://"):
    sys.exit("--post must be a https:// link (the post, or the tracking issue until the post exists)")

dest = ROOT / a.folder
dest.mkdir()
readme = (ROOT / "docs/template/README.md").read_text()
if a.title:
    readme = readme.replace("<Title: the question this folder answers>", a.title, 1)
(dest / "README.md").write_text(readme)
shutil.copy(ROOT / "docs/template/run.sh", dest / "run.sh")
(dest / "run.sh").chmod((dest / "run.sh").stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

with open(ROOT / "catalog.csv", newline="") as fh:
    rows = list(csv.DictReader(fh))
rows.append({"folder": a.folder, "area": a.area, "data": "TODO: data set and source",
             "truth": "TODO: what the results are scored against", "tools": "TODO: tools with versions",
             "post": a.post, "ci": "no"})
rows.sort(key=lambda r: r["folder"])
with open(ROOT / "catalog.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLUMNS, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)

print(f"created {a.folder}/README.md and {a.folder}/run.sh; added a catalog.csv row (fill in the TODO cells)")
print("next: write the scripts, fill in README and catalog.csv, then python3 docs/make_summaries.py and python3 ci/lint.py")
