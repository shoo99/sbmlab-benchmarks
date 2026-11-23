"""Known-answer test of single_peptide.py on synthetic search outputs (run by CI).

The same planted protein groups in the DIA-NN, MaxQuant and FragPipe formats, all passing 1%:
  targets     10 with 1 peptide, 6 with 2, 4 with 5            -> 20 target groups, 10 with one peptide (share 0.5)
  entrapment   5 with 1 peptide, 2 with 2, 1 with 3            ->  8 entrapment groups, 5 with one peptide
Left out in every format: a decoy group, a contaminant group, and a target and an entrapment group above 1%
(MaxQuant Q-value, DIA-NN Global.PG.Q.Value; FragPipe writes only what passes, so it has none). DIA-NN also gets a
precursor row above run-level q 0.01 that must not count as a peptide, a second charge state of a peptide (counted once)
and a group at 1% whose precursors all have Global.Q.Value > 0.01 (n_no_global_pep = 1).
Usage: python3 test/test_single_peptide.py   (needs pandas, as the DIA-NN reader does)
"""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
PLAN = [("T", 1)] * 10 + [("T", 2)] * 6 + [("T", 5)] * 4 + [("E", 1)] * 5 + [("E", 2)] * 2 + [("E", 3)]
EXPECT = {"n_target_groups": "20", "n_entrapment_groups": "8", "entrapment_1pep": "5", "entrapment_2pep": "2",
          "entrapment_3pep_plus": "1", "share_entrapment_1pep": "0.625", "target_1pep": "10", "share_target_1pep": "0.5",
          "n_0pep": "0"}
checks = []


def write_tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def acc(kind, i):
    return f"ENTRAP_P{i:05d}" if kind == "E" else f"P{i:05d}"


with tempfile.TemporaryDirectory() as tmp:
    t = Path(tmp)
    (t / "info.json").write_text('{"mode": "shuffle", "r_peptides": 1.0, "r_proteins": 1.0}')
    groups = [(acc(k, i), n) for i, (k, n) in enumerate(PLAN)]
    # MaxQuant
    mq = [[a, str(n), "", "", "", "", "0.001"] for a, n in groups]
    mq += [["REV__P90000", "1", "", "+", "", "", "0.001"], ["CON__P90001", "1", "", "", "+", "", "0.001"],
           ["P90002", "1", "", "", "", "", "0.02"], ["ENTRAP_P90003", "1", "", "", "", "", "0.02"]]
    write_tsv(t / "mq" / "proteinGroups.txt", ["Protein IDs", "Peptides", "Only identified by site", "Decoy",
                                               "Potential contaminant", "Reverse", "Q-value"], mq)
    # FragPipe: two experiments (each group in one of them), and the combined count
    for j in (0, 1):
        write_tsv(t / "fp" / f"exp{j}" / "protein.tsv", ["Protein", "Indistinguishable Proteins"],
                  [[f"sp|{a}|X", ""] for k, (a, n) in enumerate(groups) if k % 2 == j])
    write_tsv(t / "fp" / "combined_protein.tsv", ["Protein", "Combined Total Peptides"],
              [[f"sp|{a}|X", str(n)] for a, n in groups])
    # DIA-NN
    rows = []
    for a, n in groups:
        for p in range(n):
            rows.append([a, f"PEP{a}K{p}", "0.001", "0.001", "0.001", "0"])
        rows.append([a, f"PEP{a}K0", "0.001", "0.001", "0.001", "0"])        # a second charge state: same peptide
        rows.append([a, f"NOISE{a}K", "0.02", "0.02", "0.001", "0"])        # above run-level q: not a peptide
    rows += [["P90002", "PEPP90002K", "0.001", "0.001", "0.02", "0"], ["ENTRAP_P90003", "PEPE90003K", "0.001", "0.001", "0.02", "0"],
             ["REV__P90000", "PEPR90000K", "0.001", "0.001", "0.001", "1"], ["CON__P90001", "PEPC90001K", "0.001", "0.001", "0.001", "0"],
             ["P90004", "PEPP90004K", "0.004", "0.03", "0.001", "0"]]         # at 1%, no precursor at global 1%
    write_tsv(t / "report.tsv", ["Protein.Group", "Stripped.Sequence", "Q.Value", "Global.Q.Value", "Global.PG.Q.Value",
                                 "Decoy"], rows)
    write_tsv(t / "searches.tsv", ["search", "engine", "version", "mbr", "result", "entrapment"],
              [["mq", "maxquant", "", "off", "mq", "info.json"], ["fp", "fragpipe", "", "off", "fp", "info.json"],
               ["diann", "diann", "", "off", "report.tsv", "info.json"]])
    subprocess.run([sys.executable, str(HERE / "single_peptide.py"), str(t / "searches.tsv"), "--out", str(t / "out.csv")],
                   check=True, capture_output=True, text=True)
    got = {r["search"]: r for r in csv.DictReader(open(t / "out.csv"))}

for engine in ("mq", "fp", "diann"):
    for k, v in EXPECT.items():
        want = "21" if engine == "diann" and k == "n_target_groups" else v      # DIA-NN: + P90004
        want = "0.5238" if engine == "diann" and k == "share_target_1pep" else want
        want = "11" if engine == "diann" and k == "target_1pep" else want
        checks.append((f"{engine} {k} = {want}", got[engine][k] == want, got[engine][k]))
checks.append(("diann n_no_global_pep = 1", got["diann"]["n_no_global_pep"] == "1", got["diann"]["n_no_global_pep"]))
checks.append(("mq n_no_global_pep empty", got["mq"]["n_no_global_pep"] == "", got["mq"]["n_no_global_pep"]))
for name, ok, val in checks:
    print(f"{'ok  ' if ok else 'FAIL'} {name}" + ("" if ok else f" (got {val})"))
failed = sum(not ok for _, ok, _ in checks)
print(f"{len(checks) - failed} of {len(checks)} checks passed")
sys.exit(1 if failed else 0)
