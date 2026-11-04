"""Scores the bundled hap.py tables (happy/<caller>_<coverage>x.extended.csv) and timings.tsv into scores.csv (#14).

Call sets (hap.py Filter column):
  bcftools call -m (all calls)          bcftools rows, Filter ALL
  bcftools, QUAL >= 20                  bcftools rows, Filter PASS (soft filter LowQual on QUAL < 20)
  GATK HaplotypeCaller (unfiltered)     gatk rows, Filter ALL
  GATK HaplotypeCaller + hard filters   gatk rows, Filter PASS
  DeepVariant                           deepvariant rows, Filter PASS (its RefCall records are not calls)
Metrics per call set, variant type (SNP, INDEL) and region (all benchmark regions, and each GIAB stratum):
  recall, precision, f1 (hap.py METRIC.*), fn and fp counts, n = truth variants in the region
plus wall_seconds, cpu_seconds and peak_rss_mb of each caller (timings.tsv).
Usage: python3 score.py [happy_dir] [out_dir]   (defaults: this folder; standard library only)
"""
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
happy = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "happy"
out = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE
FOLDER = "variant-calling-giab"
COLUMNS = ["folder", "dataset", "tool", "version", "setting", "subset", "metric", "value", "n"]
SETS = [("bcftools", "ALL", "bcftools call -m (all calls)"), ("bcftools", "PASS", "bcftools, QUAL >= 20"),
        ("gatk", "ALL", "GATK HaplotypeCaller (unfiltered)"), ("gatk", "PASS", "GATK HaplotypeCaller + hard filters"),
        ("deepvariant", "PASS", "DeepVariant")]
VERSIONS = {"bcftools": "bcftools 1.25", "gatk": "GATK 4.6.2.0", "deepvariant": "DeepVariant 1.10.0"}
SUBSETS = [("*", "all benchmark regions"), ("notinalldifficultregions", "not in any difficult region"),
           ("alldifficultregions", "all difficult regions"),
           ("AllTandemRepeatsandHomopolymers_slop5", "tandem repeats and homopolymers"),
           ("segdups", "segmental duplications"), ("lowmappabilityall", "low mappability")]
SHORT = {"*": "all", "notinalldifficultregions": "easy", "alldifficultregions": "difficult",
         "AllTandemRepeatsandHomopolymers_slop5": "TR/homopol", "segdups": "segdups", "lowmappabilityall": "lowmap"}
DATASET = "GIAB HG002 chr20, NovaSeq PCR-free {cov}x, v4.2.1 benchmark (GRCh38)"


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


rows = []
covs = sorted({int(p.name.rsplit("_", 1)[1].split("x")[0]) for p in happy.glob("*_*x.extended.csv")}, reverse=True)
timings = list(csv.DictReader(open(HERE / "timings.tsv"), delimiter="\t")) if (HERE / "timings.tsv").exists() else []
for cov in covs:
    dataset = DATASET.format(cov=cov)
    print(f"== {dataset}")
    for vtype in ("SNP", "INDEL"):
        print(f"\n-- {vtype}: F1 (recall / precision) by region")
        print(f"   {'call set':<38}" + "".join(f"{SHORT[s]:>22}" for s, _ in SUBSETS))
        for caller, filt, label in SETS:
            f = happy / f"{caller}_{cov}x.extended.csv"
            if not f.exists():
                continue
            data = {(r["Type"], r["Subset"], r["Filter"]): r for r in csv.DictReader(open(f))}
            cells = []
            for subset, desc in SUBSETS:
                r = data.get((vtype, subset, filt))
                if r is None or num(r["METRIC.Recall"]) is None:
                    cells.append(f"{'-':>22}")
                    continue
                truth = int(float(r["TRUTH.TOTAL"]))
                m = {"recall": num(r["METRIC.Recall"]), "precision": num(r["METRIC.Precision"]),
                     "f1": num(r["METRIC.F1_Score"]), "fn": int(float(r["TRUTH.FN"])), "fp": int(float(r["QUERY.FP"]))}
                for metric, v in m.items():
                    if v is None:
                        continue
                    rows.append([FOLDER, dataset, label, VERSIONS[caller], f"{vtype}s", desc, metric,
                                 round(v, 5) if isinstance(v, float) else v, truth])
                cells.append(f"{m['f1']:.4f} ({m['recall']:.3f}/{m['precision'] if m['precision'] is not None else float('nan'):.3f})".rjust(22))
            print(f"   {label:<38}" + "".join(cells))
    t = [x for x in timings if x["coverage"] == str(cov)]
    if t:
        print("\n-- run time on chr20 (4 threads; peak RSS of the largest process)")
        for x in t:
            print(f"   {x['caller']:<12} wall {int(x['wall_s']):>6} s   CPU {int(x['cpu_s']):>7} s   peak RSS {x['peak_rss_mb']:>6} MB")
            for metric, col in (("wall_seconds", "wall_s"), ("cpu_seconds", "cpu_s"), ("peak_rss_mb", "peak_rss_mb")):
                rows.append([FOLDER, dataset, x["caller"], VERSIONS.get(x["caller"], x["caller"]), "chr20, 4 threads",
                             "run time", metric, int(float(x[col])), 1])
    print()
out.mkdir(parents=True, exist_ok=True)
with open(out / "scores.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(COLUMNS)
    w.writerows(rows)
print(f"wrote scores.csv ({len(rows)} rows)")
