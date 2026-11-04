"""Known-answer tests for go_ora.py (go-release-drift, issue #19): a five-term ontology, a small GAF and exact
hypergeometric and Benjamini-Hochberg values. Usage: python3 test/test_go_ora.py   (standard library only)"""
import gzip
import math
import sys
import tempfile
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from go_ora import ancestors_fn, bh, hyper_upper, log_factorials, read_gaf, read_obo  # noqa: E402

OBO = """format-version: 1.2
data-version: releases/test

[Term]
id: GO:0000001
name: root process
namespace: biological_process

[Term]
id: GO:0000002
name: child by is_a
namespace: biological_process
alt_id: GO:0000099
is_a: GO:0000001 ! root process

[Term]
id: GO:0000003
name: part of child
namespace: biological_process
relationship: part_of GO:0000002 ! child by is_a

[Term]
id: GO:0000004
name: regulates child (not propagated)
namespace: biological_process
relationship: regulates GO:0000002 ! child by is_a

[Term]
id: GO:0000005
name: obsolete thing
namespace: biological_process
is_obsolete: true

[Typedef]
id: part_of
name: part of
"""
GAF = [
    "!gaf-version: 2.2",
    "UniProtKB\tP00001\tGENE1\tinvolved_in\tGO:0000003\tPMID:1\tIDA\t\tP\t\t\tprotein\ttaxon:9606\t20200101\tX",
    "UniProtKB\tP00002\tGENE2\tNOT|involved_in\tGO:0000003\tPMID:1\tIDA\t\tP\t\t\tprotein\ttaxon:9606\t20200101\tX",
    "UniProtKB\tP00003-2\tGENE3\tinvolved_in\tGO:0000099\tGO_REF:1\tIEA\t\tP\t\t\tprotein\ttaxon:9606\t20200101\tX",
    "UniProtKB\tP00004\tGENE4\tinvolved_in\tGO:0000005\tPMID:1\tIMP\t\tP\t\t\tprotein\ttaxon:9606\t20200101\tX",
    "UniProtKB\tP00005\tGENE5\tinvolved_in\tGO:0000004\tPMID:1\tIMP\t\tP\t\t\tprotein\ttaxon:10090\t20200101\tX",
    "RNAcentral\tURS0001\tRNA1\tinvolved_in\tGO:0000004\tPMID:1\tIMP\t\tP\t\t\tRNA\ttaxon:9606\t20200101\tX",
]
failures = 0


def check(name, ok):
    global failures
    print(("ok   " if ok else "FAIL ") + name)
    failures += not ok


with tempfile.TemporaryDirectory() as d:
    obo, gaf = Path(d) / "t.obo", Path(d) / "t.gaf.gz"
    obo.write_text(OBO)
    with gzip.open(gaf, "wt") as fh:
        fh.write("\n".join(GAF) + "\n")
    terms, alt = read_obo(obo)
    check("five terms read, Typedef skipped", sorted(terms) == [f"GO:000000{i}" for i in range(1, 6)])
    check("alt_id mapped to its primary id", alt == {"GO:0000099": "GO:0000002"})
    check("obsolete flag read", terms["GO:0000005"]["obsolete"] and not terms["GO:0000001"]["obsolete"])
    anc = ancestors_fn(terms)
    check("part_of and is_a propagate", anc("GO:0000003") == {"GO:0000001", "GO:0000002", "GO:0000003"})
    check("regulates does not propagate", anc("GO:0000004") == {"GO:0000004"})
    direct = read_gaf(gaf, terms, alt, {"all": lambda ev: True, "no_IEA": lambda ev: ev != "IEA"})
    check("NOT rows, obsolete terms, other species and non-UniProtKB rows dropped",
          sorted(direct["all"]) == ["P00001", "P00003"])
    check("isoform suffix stripped and alt_id mapped", direct["all"]["P00003"] == {"GO:0000002"})
    check("IEA excluded in the no_IEA setting", sorted(direct["no_IEA"]) == ["P00001"])

lf = log_factorials(200)
check("log factorials", abs(lf[10] - math.log(math.factorial(10))) < 1e-12)


def exact_upper(k, n, K, N):
    tot = sum(Fraction(math.comb(K, i) * math.comb(N - K, n - i)) for i in range(k, min(K, n) + 1))
    return float(tot / math.comb(N, n))


worst = 0.0
for k, n, K, N in [(1, 10, 20, 100), (5, 30, 12, 150), (12, 40, 15, 200), (3, 50, 80, 200), (0, 5, 5, 50),
                   (20, 60, 25, 200), (2, 2, 2, 10)]:
    p, e = hyper_upper(k, n, K, N, lf), exact_upper(k, n, K, N) if k > 0 else 1.0
    worst = max(worst, abs(p - e) / e)
check(f"hypergeometric upper tail matches exact sums (worst relative error {worst:.1e})", worst < 1e-9)
check("impossible count gives p = 0", hyper_upper(6, 5, 10, 100, lf) == 0.0)
p = [0.01, 0.04, 0.03, 0.5, 0.2]
q = bh(p)
expected = [0.05, 0.0666666667, 0.0666666667, 0.5, 0.25]
check("Benjamini-Hochberg matches hand-computed values", all(abs(a - b) < 1e-9 for a, b in zip(q, expected)))
print("all tests passed" if not failures else f"{failures} test(s) failed")
sys.exit(1 if failures else 0)
