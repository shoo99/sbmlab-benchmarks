"""Scores the taxonomic profiles in results/ against the ZymoBIOMICS theoretical compositions (truth.tsv); issue #15.

Inputs per data set (D6300 even, D6310 log distribution) and profiler:
  results/<set>_<db>.bracken + .kreport   Kraken2 + Bracken (fraction_total_reads at species level; classified share)
  results/<set>_<tool>.profile.tsv        other profilers, converted by run.sh to: taxid, species, percent
Truth: genomic-DNA composition from the Zymo product sheets (D6300: 8 bacteria at 12%, 2 yeasts at 2%; D6310: 89.1% down
to 0.000089%). Read-based profiles estimate DNA (read) fractions, so this is the composition they should report.
Species are matched by NCBI taxid, else by name (Lactobacillus fermentum = Limosilactobacillus fermentum).
Metrics per data set and profile:
  l1_distance            sum over all species of |estimated % - true %| (0 = perfect, 200 = disjoint)
  bray_curtis            l1 / 200
  precision / recall     species called present at three thresholds: any reads reported (Bracken >= 10 reads),
                         >= 0.01 %, >= 0.1 %; a call outside the 10 species is a false positive
  false_positive_species species called (any reads) that are not in the mock
  detection_limit_percent  (D6310) the lowest true abundance among the mock species that were detected
  estimated_percent      per mock species
  classified_share       (Kraken2) share of reads classified at any rank
Usage: python3 score.py [results_dir] [out_dir]   (defaults: this folder; standard library only)
"""
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
res = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "results"
out = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE
FOLDER = "metagenome-profilers-mock"
COLUMNS = ["folder", "dataset", "tool", "version", "setting", "subset", "metric", "value", "n"]
LABEL = {"D6300": "ZymoBIOMICS D6300 (even), ERR2984773, MiSeq 2x151",
         "D6310": "ZymoBIOMICS D6310 (log distribution), ERR2935805, HiSeq 2x101"}
TOOLS = {"k2_standard_08gb": ("Kraken2 + Bracken, standard-8 (no fungi)", "Kraken2 2.1 / Bracken 3; DB 2024-09-04"),
         "k2_pluspf_08gb": ("Kraken2 + Bracken, PlusPF-8 (with fungi)", "Kraken2 2.1 / Bracken 3; DB 2024-09-04"),
         "k2_standard": ("Kraken2 + Bracken, standard (uncapped)", "Kraken2 2.1 / Bracken 3; DB 2024-09-04"),
         "metaphlan": ("MetaPhlAn 4", "see README"), "motus": ("mOTUs 3", "see README"), "sylph": ("sylph, GTDB r220", "see README")}
THRESHOLDS = [("any reads", 0.0), (">= 0.01 %", 0.01), (">= 0.1 %", 0.1)]


def norm(name):
    n = re.sub(r"[^a-z ]", "", name.lower()).strip()
    n = " ".join(n.split()[:2])
    return {"lactobacillus fermentum": "limosilactobacillus fermentum"}.get(n, n)


truth = {}
for r in csv.DictReader(open(HERE / "truth.tsv"), delimiter="\t"):
    truth.setdefault(r["dataset"], []).append(r)


def load(path):
    """[(taxid, species, percent)] summing to 100 over the species reported."""
    rows = []
    if path.suffix == ".bracken":
        for r in csv.DictReader(open(path), delimiter="\t"):
            rows.append((r["taxonomy_id"], r["name"], 100 * float(r["fraction_total_reads"])))
    else:
        for r in csv.DictReader(open(path), delimiter="\t"):
            rows.append((r["taxid"], r["species"], float(r["percent"])))
    s = sum(p for _, _, p in rows) or 1
    return [(t, n, 100 * p / s) for t, n, p in rows]


def classified_share(kreport):
    unc = root = 0
    for line in open(kreport):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 6 and f[3] == "U":
            unc = int(f[1])
        elif len(f) >= 6 and f[3] == "R" and f[5].strip() == "root":
            root = int(f[1])
    return root / (root + unc) if root + unc else None


rows = []
files = sorted(list(res.glob("*.bracken")) + list(res.glob("*.profile.tsv")))
for ds in ("D6300", "D6310"):
    t = truth[ds]
    tax_to_sp = {r["taxid"]: r["species"] for r in t}
    name_to_sp = {norm(r["species"]): r["species"] for r in t}
    true_pct = {r["species"]: float(r["gdna_percent"]) for r in t}
    ds_files = [f for f in files if f.name.startswith(ds + "_")]
    if not ds_files:
        continue
    dataset = LABEL[ds]
    print(f"== {dataset}\n")
    head = f"   {'profile':<44} {'L1':>6} {'BC':>6} {'FP spp':>6}  " + "  ".join(f"P/R {t_[0]:<9}" for t_ in THRESHOLDS)
    print(head + ("  detect.limit" if ds == "D6310" else ""))
    est_table = {}
    for f in ds_files:
        key = f.name[len(ds) + 1:].replace(".bracken", "").replace(".profile.tsv", "")
        label, version = TOOLS.get(key, (key, ""))
        prof = load(f)
        est = {sp: 0.0 for sp in true_pct}
        others = {}
        for taxid, name, p in prof:
            sp = tax_to_sp.get(taxid) or name_to_sp.get(norm(name))
            if sp:
                est[sp] += p
            else:
                others[name] = others.get(name, 0) + p
        l1 = sum(abs(est[s] - true_pct[s]) for s in true_pct) + sum(others.values())
        m = {"l1_distance": round(l1, 3), "bray_curtis": round(l1 / 200, 4),
             "false_positive_species": sum(1 for p in others.values() if p > 0)}
        pr = []
        for tname, thr in THRESHOLDS:
            tp = sum(1 for s in true_pct if est[s] > thr)
            fp = sum(1 for p in others.values() if p > thr)
            prec = tp / (tp + fp) if tp + fp else 0.0
            rec = tp / len(true_pct)
            m[f"precision ({tname})"] = round(prec, 4)
            m[f"recall ({tname})"] = round(rec, 4)
            pr.append(f"{prec:.2f}/{rec:.2f}".ljust(13))
        detected = [true_pct[s] for s in true_pct if est[s] > 0]
        if ds == "D6310":
            m["detection_limit_percent"] = min(detected) if detected else None
        kr = f.with_suffix(".kreport")
        if f.suffix == ".bracken" and kr.exists():
            cs = classified_share(kr)
            if cs is not None:
                m["classified_share"] = round(cs, 4)
        for metric, v in m.items():
            if v is None:
                continue
            base, _, thr = metric.partition(" (")
            rows.append([FOLDER, dataset, label, version, thr.rstrip(")") or "species level", "all species",
                         base, v, len(true_pct)])
        for s in true_pct:
            rows.append([FOLDER, dataset, label, version, "species level", s, "estimated_percent", round(est[s], 5), 1])
        est_table[label] = est
        dl = f"  {m['detection_limit_percent']:g} %" if ds == "D6310" and m.get("detection_limit_percent") else ""
        print(f"   {label:<44} {l1:>6.1f} {l1 / 200:>6.3f} {m['false_positive_species']:>6}  " + "  ".join(pr) + dl)
        top = sorted(others.items(), key=lambda x: -x[1])[:3]
        print(f"   {'':<44} largest false positives: " + "; ".join(f"{n} {p:.2f} %" for n, p in top))
    print(f"\n   {'species':<32} {'true %':>10}  " + "  ".join(f"{k[:22]:>22}" for k in est_table))
    for s in true_pct:
        print(f"   {s:<32} {true_pct[s]:>10g}  " + "  ".join(f"{est_table[k][s]:>22.4f}" for k in est_table))
    print()

t = list(csv.DictReader(open(HERE / "timings.tsv"), delimiter="\t")) if (HERE / "timings.tsv").exists() else []
if t:
    print("== run time (4 threads; peak RSS of the largest process; reads streamed from the .sra)")
    for x in t:
        print(f"   {x['dataset']:<6} {x['tool']:<32} wall {int(x['wall_s']):>6} s   CPU {int(x['cpu_s']):>7} s   peak RSS {x['peak_rss_mb']:>6} MB")
        for metric, col in (("wall_seconds", "wall_s"), ("cpu_seconds", "cpu_s"), ("peak_rss_mb", "peak_rss_mb")):
            rows.append([FOLDER, LABEL[x["dataset"]], x["tool"], "", "4 threads", "run time", metric, int(float(x[col])), 1])
out.mkdir(parents=True, exist_ok=True)
with open(out / "scores.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(COLUMNS)
    w.writerows(rows)
print(f"\nwrote scores.csv ({len(rows)} rows)")
