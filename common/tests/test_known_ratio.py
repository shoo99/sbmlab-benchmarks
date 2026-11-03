"""Tests for common/known_ratio.py: the same toy matrix as test_known_ratio.R, then re-creating the bundled
diann-hye/diann_scores.csv and diann-2.7-hye/diann_scores_270.csv. Run: python3 common/tests/test_known_ratio.py
"""
import csv
import math
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "common"))
import known_ratio as kr  # noqa: E402

fails = 0


def check(desc, ok):
    global fails
    print(("ok   " if ok else "FAIL ") + desc)
    fails += not ok


def close(a, b):
    return all((math.isnan(x) and math.isnan(y)) or abs(x - y) < 1e-9 for x, y in zip(a, b)) and len(a) == len(b)


L = np.array([[10, 10, 10, 10, 10, 10], [10, 10, 10, 10.4, 10.4, 10.4], [12, 12, 12, 11, 11, 11],
              [8, 8, np.nan, 10.6, 10.6, 10.6], [8, np.nan, np.nan, 10, 10, 10], [5] * 6, [5] * 6,
              [10, 10, 10, 11, 11, 11]])
m = np.nan_to_num(2.0 ** L, nan=0.0)
nm = ["H1_HUMAN", "H2_HUMAN", "Y1_YEAST", "E1_ECOLI", "E2_ECOLI", "X_HUMAN;Y_YEAST", "Z_HUMANX", "Y2_YEAST"]

sp = kr.species_from_names(nm)
check("species_from_names", sp == ["Human", "Human", "Yeast", "Ecoli", "Ecoli", None, None, "Yeast"])
check("species_from_taxa", kr.species_from_taxa(["OX=9606", "OX=9606 OX=559292", "OX=10090", "OX=10090;OX=83333", ""])
      == ["Human", None, None, "Ecoli", None])

ok, fc = kr.log2_ratio(m, [0, 1, 2], [3, 4, 5])
check("log2_ratio ok (2 of 3 per group)", ok.tolist() == [True, True, True, True, False, True, True, True])
check("log2_ratio fc", close(fc.tolist(), [0, 0.4, -1, 2.6, 2, 0, 0, 1]))

g = kr.group_log2_means(m, {"A": [0, 1, 2], "B": [3, 4, 5]})
check("group_log2_means", close(g["A"][3:5].tolist(), [8, float("nan")]) and close([g["B"][1]], [10.4]))

rows, shift = kr.score_groups(fc, ok, sp, kr.HYE_TRUTH, centre="Human")
check("score_groups n", [r["n"] for r in rows] == [2, 2, 1])
check("score_groups median", close([r["median"] for r in rows], [0.2, 0, 2.6]))
check("score_groups IQR", close([r["IQR"] for r in rows], [0.2, 1, 0]))
check("score_groups within", close([r["within"] for r in rows], [1, 0.5, 0]))
check("score_groups within_centred", close([r["within_centred"] for r in rows], [1, 0.5, 1]))
check("score_groups centre median", close([shift], [0.2]))

o = kr.call_outcome([True, False, True, True, True, True, False, True], fc, sp, kr.HYE_TRUTH)
check("call_outcome", o == ["false", "not_called", "right", "right", "right", "no_truth", "not_called", "wrong_direction"])
cs = kr.call_summary(o)
check("call_summary", (cs["calls"], cs["right"], cs["wrong_direction"], cs["false"]) == (6, 3, 1, 1) and close([cs["fdp"]], [0.4]))

with tempfile.NamedTemporaryFile("w+", suffix=".csv", delete=False) as fh:
    path = fh.name
kr.write_scores([dict(folder="f", dataset="d", tool="t", version="1.0", setting="default", subset="Human",
                      metric="within_0.5", value=1, n=2)], path)
with open(path, newline="") as fh:
    check("write_scores header", next(csv.reader(fh)) == kr.SCORES_COLUMNS)


def read_matrix(path):
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        rows = list(rd)
    cols = [i for i, h in enumerate(header) if "Condition_" in h]
    a = [i for i in cols if "Condition_A" in header[i]]
    b = [i for i in cols if "Condition_B" in header[i]]
    names = [r[header.index("Protein.Names")] for r in rows]
    vals = np.array([[float(r[i]) if r[i] not in ("", "NA") else np.nan for i in a + b] for r in rows])
    return names, vals, len(a)


def recreate(folder, files, labels, ref):
    got = []
    for f, label in zip(files, labels):
        names, vals, na = read_matrix(ROOT / folder / f)
        ok, fc = kr.log2_ratio(vals, range(na), range(na, vals.shape[1]))
        rows, _ = kr.score_groups(fc, ok, kr.species_from_names(names), kr.HYE_TRUTH, centre="Human")
        got += [[label, r["group"], r["n"], r["median"], r["IQR"], r["within"], r["within_centred"]] for r in rows]
    with open(ROOT / folder / ref, newline="") as fh:
        bundled = [[r[0], r[1], int(r[2])] + [float(x) for x in r[3:]] for r in list(csv.reader(fh))[1:]]
    check(f"re-creates {folder}/{ref}", got == bundled)


recreate("diann-hye", ["report.pg_matrix.tsv", "pg_mbr.tsv", "pg_matrix_nombr_usequant.tsv"],
         ["pg_matrix (MBR, DIA-NN default filter)", "MBR pass, Q<=0.01 & PG.Q<=0.01", "no MBR (first-pass .quant reused)"],
         "diann_scores.csv")
recreate("diann-2.7-hye", ["report.pg_matrix.tsv", "pg_mbr.tsv"],
         ["pg_matrix (MBR, DIA-NN default filter)", "MBR pass, Q<=0.01 & PG.Q<=0.01"], "diann_scores_270.csv")

if fails:
    print(f"{fails} test(s) failed")
    sys.exit(1)
print("all Python tests passed")
