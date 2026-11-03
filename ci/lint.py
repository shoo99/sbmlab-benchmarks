"""Repository checks that need only the Python standard library.

Run by .github/workflows/checks.yml (job "static"); run locally with `python3 ci/lint.py`.

1. catalog.csv has one row per analysis folder and no other rows, with every column filled in.
2. Folders marked ci=reproduce or ci=smoke in catalog.csv are run by ci/reproduce.sh, and ci=no folders are not.
3. Every tracked .json file is strict JSON (no NaN or Infinity, which many parsers reject).
4. Every tracked .py file parses, and every tracked .sh file passes `bash -n`.
5. No script writes to an absolute path. Outputs belong next to the script; inputs may live anywhere.
   ALLOW lists the exceptions and why.
6. Every scores.csv has the standard header (see CONTRIBUTING.md).
"""
import ast
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NOT_ANALYSES = {"docs", "ci"}                     # top-level folders that are not analysis folders
COLUMNS = ["folder", "area", "data", "truth", "tools", "post", "ci"]
CI_VALUES = {"reproduce", "smoke", "no"}
SCORES_HEADER = ["folder", "dataset", "tool", "version", "setting", "subset", "metric", "value", "n"]
ALLOW = {"string-hubs/figure.py": "writes the post's figure into the blog repository"}

R_WRITE = re.compile(r"\b(fwrite|write\.csv|write\.table|writeLines|saveRDS|ggsave|png|pdf|svg|sink)\s*\(")
PY_WRITE = re.compile(r"\b(to_csv|to_parquet|to_excel|to_json|savefig)\s*\(|\bopen\([^)]*,\s*[\"'][wax]")
ABSOLUTE = re.compile(r"[\"'](/home/|/mnt/|/Users/|/tmp/|~/|[A-Za-z]:\\\\)")

errors = []


def git_files():
    out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
    return [line for line in out.splitlines() if line]


def check_catalog(files):
    folders = {f.split("/")[0] for f in files if "/" in f and not f.startswith(".")} - NOT_ANALYSES
    with open(ROOT / "catalog.csv", newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    if reader.fieldnames != COLUMNS:
        errors.append(f"catalog.csv: header must be {','.join(COLUMNS)}")
        return
    names = [r["folder"] for r in rows]
    errors.extend(f"catalog.csv: folder {d}/ has no row" for d in sorted(folders - set(names)))
    errors.extend(f"catalog.csv: row '{d}' has no folder" for d in sorted(set(names) - folders))
    errors.extend(f"catalog.csv: '{d}' is listed more than once" for d in sorted({n for n in names if names.count(n) > 1}))
    reproduce = (ROOT / "ci" / "reproduce.sh").read_text()
    for r in rows:
        empty = [c for c in COLUMNS if not (r[c] or "").strip()]
        if empty:
            errors.append(f"catalog.csv: {r['folder']}: empty {', '.join(empty)}")
        if r["ci"] not in CI_VALUES:
            errors.append(f"catalog.csv: {r['folder']}: ci must be one of {', '.join(sorted(CI_VALUES))}")
        if not r["post"].startswith("https://"):
            errors.append(f"catalog.csv: {r['folder']}: post must be a https:// link")
        run = re.search(rf"^\s*(check_out|run_ok|same_file)\s+{re.escape(r['folder'])}[/\s]", reproduce, re.M)
        if r["ci"] in ("reproduce", "smoke") and not run:
            errors.append(f"catalog.csv: {r['folder']} has ci={r['ci']} but ci/reproduce.sh does not run it")
        if r["ci"] == "no" and run:
            errors.append(f"catalog.csv: {r['folder']} has ci=no but ci/reproduce.sh runs it")


def no_nan(token):
    raise ValueError(f"non-standard JSON constant {token}")


def check_files(files):
    for f in files:
        path = ROOT / f
        if f.endswith(".json"):
            try:
                json.loads(path.read_text(), parse_constant=no_nan)
            except ValueError as e:
                errors.append(f"{f}: not strict JSON ({e})")
        elif f.endswith(".py"):
            try:
                ast.parse(path.read_text(), filename=f)
            except SyntaxError as e:
                errors.append(f"{f}:{e.lineno}: syntax error ({e.msg})")
        elif f.endswith(".sh"):
            r = subprocess.run(["bash", "-n", str(path)], capture_output=True, text=True)
            if r.returncode:
                errors.append(f"{f}: bash -n failed: {r.stderr.strip()}")
        if f.endswith((".R", ".py")) and f not in ALLOW:
            writes = R_WRITE if f.endswith(".R") else PY_WRITE
            for n, line in enumerate(path.read_text().splitlines(), 1):
                if writes.search(line) and ABSOLUTE.search(line):
                    errors.append(f"{f}:{n}: writes to an absolute path; write next to the script or take the path as an argument")
        if Path(f).name == "scores.csv":
            with open(path, newline="") as fh:
                header = next(csv.reader(fh), [])
            if header != SCORES_HEADER:
                errors.append(f"{f}: header must be {','.join(SCORES_HEADER)}")


files = git_files()
check_catalog(files)
check_files(files)
if errors:
    print("\n".join(errors))
    print(f"\n{len(errors)} problem(s)")
    sys.exit(1)
print(f"ok: catalog.csv, {len(files)} tracked files")
