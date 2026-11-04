"""Scores the doublet detectors from cell_scores.tsv.gz and datasets.tsv (written by collect.py).

Only droplets whose truth is known (doublet or singlet) are scored; "unknown" droplets are ranked but not counted.
  auprc                  average precision of the score (ties grouped); the baseline is the doublet share
  doublet_share          that baseline: known doublets / known droplets
  recall_at_expected     the top k droplets by score are called, k = expected rate x all droplets;
  precision_at_expected  recall and precision of those calls among the known droplets
  recall_own_call        the tool's own doublet calls, where it makes them
  precision_own_call
Precision counts calls on invisible doublets (same donor or same species) as false, so it is a lower bound;
datasets.tsv gives the share of doublets each truth can see.
Usage: python3 score.py [bundle_dir] [out_dir]   (both default to this folder; standard library only)
Writes scores.csv (standard columns) into out_dir and prints a table per data set.
"""
import csv
import gzip
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
bundle = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
out = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE
FOLDER = "scrna-doublets"
COLUMNS = ["folder", "dataset", "tool", "version", "setting", "subset", "metric", "value", "n"]
LABEL = {"kang_ctrl": "Kang 2018 control PBMC (demuxlet truth)", "kang_stim": "Kang 2018 IFN-b PBMC (demuxlet truth)",
         "hgmm_12k": "10x hgmm_12k human + mouse (species truth)"}

info = {r["dataset"]: r for r in csv.DictReader(open(bundle / "datasets.tsv"), delimiter="\t")}
with gzip.open(bundle / "cell_scores.tsv.gz", "rt") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    tools = [c[:-6] for c in reader.fieldnames if c.endswith(" score")]
    data = {}
    for r in reader:
        data.setdefault(r["dataset"], []).append(r)


def average_precision(scores, labels):
    """Area under the precision-recall curve as a step function over distinct thresholds (ties form one step)."""
    pairs = sorted(zip(scores, labels), key=lambda x: -x[0])
    pos = sum(labels)
    ap, tp, fp, i = 0.0, 0, 0, 0
    while i < len(pairs):
        j = i
        while j < len(pairs) and pairs[j][0] == pairs[i][0]:
            tp += pairs[j][1]
            fp += 1 - pairs[j][1]
            j += 1
        new_tp = sum(p[1] for p in pairs[i:j])
        if new_tp:
            ap += (new_tp / pos) * (tp / (tp + fp))
        i = j
    return ap


def versions(ds):
    out = {}
    for part in info[ds]["tools"].split("; "):
        tool, _, rest = part.partition(" = ")
        out[tool] = rest
    return out


rows = []
for ds, cells in data.items():
    ver = versions(ds)
    n_all = len(cells)
    k = round(float(info[ds]["expected_rate"]) * n_all)
    known = [c for c in cells if c["truth"] in ("doublet", "singlet")]
    pos = sum(c["truth"] == "doublet" for c in known)
    dataset = f"{LABEL.get(ds, ds)}: {n_all} droplets"
    print(f"\n== {dataset}; {len(known)} with known truth ({pos} doublets, share {pos / len(known):.3f}); "
          f"expected rate {float(info[ds]['expected_rate']):.3f} -> top {k} called; truth sees "
          f"{100 * float(info[ds]['visible_share']):.0f}% of doublets")
    print(f"   {'tool':<30} {'auprc':>6} {'rec@exp':>8} {'prec@exp':>9} {'rec_own':>8} {'prec_own':>9}")
    for t in tools:
        if not cells[0][f"{t} score"]:
            continue
        order = sorted(range(n_all), key=lambda i: (-float(cells[i][f"{t} score"]), i))
        called = set(order[:k])
        idx_known = [i for i in range(n_all) if cells[i]["truth"] in ("doublet", "singlet")]
        lab = [1 if cells[i]["truth"] == "doublet" else 0 for i in idx_known]
        sc = [float(cells[i][f"{t} score"]) for i in idx_known]
        m = {"auprc": average_precision(sc, lab), "doublet_share": pos / len(known)}
        tp = sum(1 for i, y in zip(idx_known, lab) if y and i in called)
        ck = sum(1 for i in idx_known if i in called)
        m["recall_at_expected"] = tp / pos
        m["precision_at_expected"] = tp / ck if ck else 0.0
        own = [cells[i][f"{t} call"] for i in idx_known]
        if all(o in ("0", "1") for o in own):
            tpo = sum(1 for o, y in zip(own, lab) if o == "1" and y)
            co = sum(1 for o in own if o == "1")
            m["recall_own_call"] = tpo / pos
            m["precision_own_call"] = tpo / co if co else 0.0
        version, _, setting = ver.get(t, "").partition(" (")
        for metric, v in m.items():
            rows.append([FOLDER, dataset, t, version, setting.rstrip(")") or "default", "droplets with known truth",
                         metric, round(v, 4), len(known)])
        fmt = lambda key: f"{m[key]:.3f}" if key in m else "-"
        print(f"   {t:<30} {fmt('auprc'):>6} {fmt('recall_at_expected'):>8} {fmt('precision_at_expected'):>9} "
              f"{fmt('recall_own_call'):>8} {fmt('precision_own_call'):>9}")
out.mkdir(parents=True, exist_ok=True)
with open(out / "scores.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(COLUMNS)
    w.writerows(rows)
print(f"\nwrote scores.csv ({len(rows)} rows)")
