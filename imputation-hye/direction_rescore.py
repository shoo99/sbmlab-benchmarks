"""Direction-aware re-score of the bundled imputation results (stdlib only).

impute_bench.R before the 2026-11-02 fix counted every non-human call as a "changed protein found", including
yeast called up and E. coli called down, and left those wrong-direction calls out of the FDP. The bundled
impute_<tag>.out / .csv.gz come from that version. The .csv.gz keeps the centred log2FC (lfc_c = logFC - human
median) and adj.P, but not the human median itself, so this script recovers it: for each method it finds the
offsets h for which
    call = adj.P < 0.05 & |lfc_c + h| >= 0.5
reproduces the "changed found" and "human FP" counts printed in impute_<tag>.out, then counts the calls by
direction (logFC = lfc_c + h; yeast must go down, E. coli up), and checks that the corrected counts are the same
for every matching offset. (Cross-check: the offset found for 'none', MBR on, is +0.144, the human median of
MaxQuant LFQ in maxquant-hye/quant_scores.csv.)

Writes impute_<tag>_summary.csv in the format the fixed impute_bench.R writes, and prints the table
(direction_rescore.out). If impute_<tag>.out already comes from the fixed script, there is nothing to re-score.

Usage: cd imputation-hye && python3 direction_rescore.py > direction_rescore.out
"""
import csv, gzip, re, sys
from collections import defaultdict

ROW = re.compile(r"^(\S+)\s+\|\s+(\d+)\s+\|.*\|\s+(\d+) of\s+(\d+)\s+\|\s+(\d+)\s+\|\s+([\d.]+)%")
FILTERS = {"BOTH": ">=2 of 3 in BOTH groups", "EITHER": ">=2 of 3 in EITHER group"}


def printed(tag):
    """[(filter, method, testable, found, n_changed, human_fp)] as printed by the original impute_bench.R."""
    out, filt = [], None
    for line in open(f"impute_{tag}.out"):
        if "wrong dir" in line:
            sys.exit(f"impute_{tag}.out comes from the fixed impute_bench.R (it has a 'wrong dir' column); "
                     f"its impute_{tag}_summary.csv is already direction-aware.")
        if line.startswith("===="):
            filt = "BOTH" if "BOTH" in line else "EITHER"
        m = ROW.match(line)
        if m and filt:
            out.append((filt, m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4)), int(m.group(5))))
    return out


def rows(tag):
    d = defaultdict(list)
    with gzip.open(f"impute_{tag}.csv.gz", "rt") as fh:
        for r in csv.DictReader(fh):
            if r["adjp"] in ("", "NA"):
                continue
            f = "BOTH" if "BOTH" in r["filter"] else "EITHER"
            d[(f, r["method"])].append((r["species"], float(r["lfc_c"]), float(r["adjp"])))
    return d


def counts(rs, h):
    """(non-human calls, human calls, right-direction calls, wrong-direction calls) at offset h."""
    found_any = fp = right = wrong = 0
    for s, lc, p in rs:
        lfc = lc + h
        if not (p < 0.05 and abs(lfc) >= 0.5):
            continue
        if s == "Human":
            fp += 1
            continue
        found_any += 1
        ok = (s == "Yeast" and lfc < 0) or (s == "Ecoli" and lfc > 0)
        right += ok
        wrong += not ok
    return found_any, fp, right, wrong


def offsets(rs, lo=-1.5, hi=1.5):
    """One candidate offset per interval between the breakpoints where some |lfc_c + h| crosses 0.5."""
    bp = sorted({b for _, lc, p in rs if p < 0.05 for b in (0.5 - lc, -0.5 - lc) if lo < b < hi} | {lo, hi})
    return [(a + b) / 2 for a, b in zip(bp, bp[1:])]


print("data   | filter | method  | printed found / FP / FDP | offsets        | correct dir   | wrong dir | human FP | FDP (human + wrong dir)")
for tag in ("MBRon", "MBRoff"):
    pr = printed(tag)
    d, summary = rows(tag), []
    for filt, meth, testable, found, n_changed, fp_pr in pr:
        rs = d[(filt, meth)]
        hits = [(h, c) for h in offsets(rs) for c in [counts(rs, h)] if c[0] == found and c[1] == fp_pr]
        res = {c for _, c in hits}
        if len(res) != 1:
            sys.exit(f"{tag} {filt} {meth}: {'no offset reproduces' if not res else 'offsets disagree on'} the printed counts")
        _, fp, right, wrong = res.pop()
        hs = [h for h, _ in hits]
        calls = found + fp_pr
        fdp = 100 * (fp + wrong) / max(calls, 1)
        print(f"{tag:6} | {filt:6} | {meth:7} | {found:5} / {fp_pr:3} / {100 * fp_pr / max(calls, 1):4.1f}%      | "
              f"{min(hs):+.3f}..{max(hs):+.3f} | {right:5} of {n_changed:4} | {wrong:9} | {fp:8} | {fdp:4.1f}%")
        summary.append(dict(data=tag, filter=FILTERS[filt], method=meth, testable=testable, correct=right, of=n_changed,
                            human_fp=fp, wrong_dir=wrong, calls=calls, fdp=round(fdp, 1)))
    with open(f"impute_{tag}_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]))
        w.writeheader()
        w.writerows(summary)
