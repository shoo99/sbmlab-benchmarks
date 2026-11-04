"""Counts target and entrapment identifications in the search outputs listed in a manifest (issue #5).

Manifest (tab-separated, paths relative to the manifest; see searches.example.tsv):
  search       a name, e.g. "DIA-NN 2.3 library-free, MBR on"
  engine       diann | maxquant | fragpipe
  version      the engine version, e.g. "DIA-NN 2.3.0"
  mbr          on | off
  result       DIA-NN: report.tsv or report.parquet; MaxQuant: the combined/txt folder; FragPipe: the output folder
  baseline     optional: the same search against the target database alone (for the targets lost)
  entrapment   entrapment_info.json of the database used (from make_entrapment.py)
Levels and q-values
  diann     precursor (Precursor.Id; Global.Q.Value, else the lowest Q.Value over runs), peptide (Stripped.Sequence;
            lowest q of its precursors), protein (Protein.Group; Global.PG.Q.Value, else the lowest PG.Q.Value)
  maxquant  precursor (evidence.txt: modified sequence + charge), peptide (peptides.txt), protein (proteinGroups.txt);
            MaxQuant writes only identifications passing its FDR, so the threshold is "reported" (Q-value is used for
            protein groups when the column is there)
  fragpipe  precursor (psm.tsv: modified peptide + charge), peptide (peptide.tsv), protein (protein.tsv), all folders
            below result; Philosopher filters at its FDR, so the threshold is "reported"
An identification is entrapment when every protein it maps to is an entrapment protein (accession contains ENTRAP_),
decoy (rev_ / REV__ / DECOY) or contaminant (CON__ / contam_) when all its proteins are; otherwise target. Decoys and
contaminants are not counted.
Output: counts.csv next to this script (or --out): search, engine, version, mbr, level, threshold, n_target,
n_entrapment, n_baseline_target, r. Thresholds: 0.001, 0.005, 0.01, 0.02, 0.05 where q-values exist, else "reported".
Usage: python3 extract.py <searches.tsv> [--out counts.csv]
Standard library only; DIA-NN .parquet reports need pandas with pyarrow.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

GRID = [0.001, 0.005, 0.01, 0.02, 0.05]
csv.field_size_limit(sys.maxsize)


def protein_class(accessions):
    acc = [a for a in accessions if a]
    if not acc:
        return "other"
    if all(a.startswith(("rev_", "REV__", "DECOY", "decoy_")) for a in acc):
        return "decoy"
    if all(a.startswith(("CON__", "contam_", "Cont_")) for a in acc):
        return "contaminant"
    real = [a for a in acc if not a.startswith(("rev_", "REV__", "DECOY", "decoy_", "CON__", "contam_", "Cont_"))]
    if not real:
        return "other"
    return "entrapment" if all("ENTRAP_" in a for a in real) else "target"


def split_proteins(s):
    return [p.strip() for p in s.replace(",", ";").split(";") if p.strip()]


def rows_tsv(path):
    with open(path, newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def best(d, key, q, prot):
    old = d.get(key)
    if old is None or (q is not None and (old[0] is None or q < old[0])):
        d[key] = (q, prot)


def read_diann(path):
    path = Path(path)
    if path.suffix == ".parquet":
        import pandas as pd                                         # only needed for parquet reports
        rows = (dict(zip(r.index, ("" if v is None else str(v) for v in r.values))) for _, r in pd.read_parquet(path).iterrows())
    else:
        rows = rows_tsv(path)
    prec, pep, prot = {}, {}, {}
    for r in rows:
        if r.get("Decoy", "0") not in ("0", "", "0.0", "False"):
            continue
        pids = split_proteins(r.get("Protein.Ids", ""))
        q = float(r["Global.Q.Value"]) if r.get("Global.Q.Value") else float(r["Q.Value"])
        best(prec, r["Precursor.Id"], q, pids)
        best(pep, r["Stripped.Sequence"], q, pids)
        pq = r.get("Global.PG.Q.Value") or r.get("PG.Q.Value")
        best(prot, r["Protein.Group"], float(pq) if pq else None, split_proteins(r["Protein.Group"]))
    return {"precursor": prec, "peptide": pep, "protein": prot}


def read_maxquant(folder):
    folder = Path(folder)
    # MaxQuant 2.8 names the decoy column "Decoy"; older versions "Reverse"
    keep = lambda r: "+" not in (r.get("Reverse", ""), r.get("Decoy", ""), r.get("Potential contaminant", ""))
    prec, pep, prot = {}, {}, {}
    for r in rows_tsv(folder / "evidence.txt"):
        if keep(r) and r.get("Type", "") != "MULTI-MATCH-NOVERIFY":
            best(prec, (r["Modified sequence"], r["Charge"]), None, split_proteins(r["Proteins"]))
    for r in rows_tsv(folder / "peptides.txt"):
        if keep(r):
            best(pep, r["Sequence"], None, split_proteins(r["Proteins"]))
    for r in rows_tsv(folder / "proteinGroups.txt"):
        if keep(r) and r.get("Only identified by site", "") != "+":
            q = r.get("Q-value")
            best(prot, r["Protein IDs"], float(q) if q not in (None, "") else None, split_proteins(r["Protein IDs"]))
    return {"precursor": prec, "peptide": pep, "protein": prot}


def read_fragpipe(folder):
    folder = Path(folder)
    prec, pep, prot = {}, {}, {}
    mapped = lambda r, col: [r[col]] + split_proteins(r.get("Mapped Proteins", ""))
    for f in sorted(folder.rglob("psm.tsv")):
        for r in rows_tsv(f):
            best(prec, (r.get("Modified Peptide") or r["Peptide"], r["Charge"]), None, mapped(r, "Protein"))
    for f in sorted(folder.rglob("peptide.tsv")):
        for r in rows_tsv(f):
            best(pep, r["Peptide"], None, mapped(r, "Protein"))
    for f in sorted(folder.rglob("protein.tsv")):
        for r in rows_tsv(f):
            best(prot, r["Protein"], None, [r["Protein"]] + split_proteins(r.get("Indistinguishable Proteins", "")))
    return {"precursor": prec, "peptide": pep, "protein": prot}


READERS = {"diann": read_diann, "maxquant": read_maxquant, "fragpipe": read_fragpipe}


def count(levels):
    """{level: {threshold: (n_target, n_entrapment)}}"""
    out = {}
    for level, ids in levels.items():
        qs = [v[0] for v in ids.values()]
        has_q = bool(ids) and all(q is not None for q in qs)
        thresholds = GRID if has_q else ["reported"]
        out[level] = {}
        for t in thresholds:
            nt = ne = 0
            for q, prots in ids.values():
                if t != "reported" and q > t:
                    continue
                c = protein_class(prots)
                nt += c == "target"
                ne += c == "entrapment"
            out[level][t] = (nt, ne)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("manifest")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "counts.csv"))
    a = ap.parse_args()
    base = Path(a.manifest).resolve().parent
    out_rows = []
    for s in rows_tsv(a.manifest):
        info = json.loads((base / s["entrapment"]).read_text())
        res = count(READERS[s["engine"]](base / s["result"]))
        bl = count(READERS[s["engine"]](base / s["baseline"])) if s.get("baseline") else None
        for level, by_t in res.items():
            r = info["r_proteins"] if level == "protein" else info["r_peptides"]
            for t, (nt, ne) in by_t.items():
                nb = bl[level].get(t, (None, None))[0] if bl else None
                out_rows.append([s["search"], s["engine"], s.get("version", ""), s["mbr"], level, t, nt, ne,
                                 "" if nb is None else nb, r])
                print(f"{s['search']:<45} {level:<9} q<={t!s:<8} target {nt:>7}  entrapment {ne:>6}"
                      + ("" if nb is None else f"  (target-only search: {nb})"))
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["search", "engine", "version", "mbr", "level", "threshold", "n_target", "n_entrapment",
                    "n_baseline_target", "r"])
        w.writerows(out_rows)
    print(f"wrote {Path(a.out).name}: {len(out_rows)} rows")


if __name__ == "__main__":
    main()
