"""Protein groups by number of peptides at the engines' 1% FDR: how many false (entrapment) groups rest on one peptide,
and what a two-peptide rule would cost in target groups (issue #5).

Protein groups and their classes are those of extract.py at the nominal 1%: DIA-NN groups with Global.PG.Q.Value
(else the lowest PG.Q.Value) <= 0.01, MaxQuant proteinGroups.txt with Q-value <= 0.01 (decoys, contaminants and
"only identified by site" out), FragPipe protein.tsv as reported (the union over experiments). A group is entrapment
when every protein in it is ENTRAP_ (decoys and contaminants ignored), otherwise target; groups that are all decoy or
all contaminant are not counted.
Number of peptides of a group = distinct peptide sequences (no modifications, no charge) that pass the engine's 1% FDR
and belong to the group:
  maxquant  the "Peptides" column of proteinGroups.txt (peptide sequences of all proteins in the group, shared ones
            included)
  fragpipe  "Combined Total Peptides" of combined_protein.tsv (distinct peptides over all experiments that map to the
            group's protein, shared ones included); if that file is missing, the same count from every peptide.tsv
            ("Protein" and "Mapped Proteins")
  diann     distinct Stripped.Sequence of the report's non-decoy rows with that Protein.Group and a run-level
            precursor Q.Value <= 0.01 (the filter DIA-NN applies to the report at --qvalue 0.01; a precursor identified
            at 1% in at least one run, as MaxQuant and FragPipe count a peptide). DIA-NN assigns each precursor to one
            protein group, so a peptide shared between groups counts only for the group it was assigned to. The
            experiment-wide Global.Q.Value is not used: DIA-NN passes some groups at Global.PG.Q.Value <= 0.01 whose
            precursors all have Global.Q.Value > 0.01, and those groups would have no peptide (n_no_global_pep).
Output columns (one row per search):
  n_target_groups, n_entrapment_groups   as in counts.csv, level protein, at 1%
  entrapment_1pep, entrapment_2pep, entrapment_3pep_plus   entrapment groups by number of peptides
  share_entrapment_1pep                  entrapment_1pep / n_entrapment_groups
  target_1pep                            target groups with one peptide
  share_target_1pep                      target_1pep / n_target_groups: the share of target groups that a rule
                                         "at least two peptides" would remove
  n_0pep                                 groups with no peptide passing the filters above (0 expected)
  n_no_global_pep                        DIA-NN only: groups at 1% with no precursor at Global.Q.Value <= 0.01
Usage: python3 single_peptide.py <searches.tsv> [<searches.tsv> ...] [--out single_peptide.csv]
       (manifests as for extract.py; the search outputs are not in the repository). DIA-NN .parquet needs pandas.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract import protein_class, rows_tsv, split_proteins  # noqa: E402

Q = 0.01


def groups_maxquant(folder):
    out = {}
    for r in rows_tsv(Path(folder) / "proteinGroups.txt"):
        if "+" in (r.get("Reverse", ""), r.get("Decoy", ""), r.get("Potential contaminant", "")) \
                or r.get("Only identified by site", "") == "+":
            continue
        if r.get("Q-value") not in (None, "") and float(r["Q-value"]) > Q:
            continue
        out[r["Protein IDs"]] = (split_proteins(r["Protein IDs"]), int(r["Peptides"]))
    return out, ""


def groups_fragpipe(folder):
    folder = Path(folder)
    prots = {}
    for f in sorted(folder.rglob("protein.tsv")):
        for r in rows_tsv(f):
            prots.setdefault(r["Protein"], [r["Protein"]] + split_proteins(r.get("Indistinguishable Proteins", "")))
    combined = folder / "combined_protein.tsv"
    if combined.exists():
        npep = {r["Protein"]: int(r["Combined Total Peptides"]) for r in rows_tsv(combined)}
    else:
        peps = {}
        for f in sorted(folder.rglob("peptide.tsv")):
            for r in rows_tsv(f):
                for p in [r["Protein"]] + split_proteins(r.get("Mapped Proteins", "")):
                    peps.setdefault(p, set()).add(r["Peptide"])
        npep = {p: len(s) for p, s in peps.items()}
    missing = [k for k in prots if k not in npep]
    if missing:
        sys.exit(f"{folder}: {len(missing)} protein groups without a peptide count, e.g. {missing[0]}")
    return {k: (v, npep[k]) for k, v in prots.items()}, ""


def groups_diann(path):
    import pandas as pd
    path = Path(path)
    cols = ["Protein.Group", "Stripped.Sequence", "Q.Value", "Global.Q.Value", "PG.Q.Value", "Global.PG.Q.Value", "Decoy"]
    d = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, sep="\t")
    d = d[[c for c in cols if c in d.columns]]
    if "Decoy" in d.columns:
        d = d[d["Decoy"].astype(str).isin(["0", "0.0", "False", ""])]
    pq = d["Global.PG.Q.Value"] if "Global.PG.Q.Value" in d.columns else d["PG.Q.Value"]
    pgq = pq.groupby(d["Protein.Group"]).min()
    groups = pgq[pgq <= Q].index
    npep = d[d["Q.Value"] <= Q].groupby("Protein.Group")["Stripped.Sequence"].nunique()
    if "Global.Q.Value" in d.columns:
        with_global = set(d.loc[d["Global.Q.Value"] <= Q, "Protein.Group"])
        no_global = sum(g not in with_global for g in groups)
    else:
        no_global = ""
    return {g: (split_proteins(g), int(npep.get(g, 0))) for g in groups}, no_global


READ = {"diann": groups_diann, "maxquant": groups_maxquant, "fragpipe": groups_fragpipe}
COLUMNS = ["search", "engine", "version", "mbr", "entrapment", "n_target_groups", "n_entrapment_groups",
           "entrapment_1pep", "entrapment_2pep", "entrapment_3pep_plus", "share_entrapment_1pep", "target_1pep",
           "share_target_1pep", "n_0pep", "n_no_global_pep"]


def summarise(groups, no_global=""):
    t = {1: 0, 2: 0, 3: 0, 0: 0}
    e = {1: 0, 2: 0, 3: 0, 0: 0}
    for prots, n in groups.values():
        c = protein_class(prots)
        if c not in ("target", "entrapment"):
            continue
        (t if c == "target" else e)[min(n, 3)] += 1
    nt, ne = sum(t.values()), sum(e.values())
    return {"n_target_groups": nt, "n_entrapment_groups": ne, "entrapment_1pep": e[1], "entrapment_2pep": e[2],
            "entrapment_3pep_plus": e[3], "share_entrapment_1pep": round(e[1] / ne, 4) if ne else 0.0,
            "target_1pep": t[1], "share_target_1pep": round(t[1] / nt, 4) if nt else 0.0, "n_0pep": t[0] + e[0],
            "n_no_global_pep": no_global}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("manifest", nargs="+")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "single_peptide.csv"))
    a = ap.parse_args()
    rows = []
    for m in a.manifest:
        base = Path(m).resolve().parent
        for s in rows_tsv(m):
            info = json.loads((base / s["entrapment"]).read_text())
            c = summarise(*READ[s["engine"]](base / s["result"]))
            rows.append({"search": s["search"], "engine": s["engine"], "version": s.get("version", ""), "mbr": s["mbr"],
                         "entrapment": info.get("mode", ""), **c})
            print(f"{s['search'][:52]:<52} T {c['n_target_groups']:>5} E {c['n_entrapment_groups']:>3}  E by peptides "
                  f"1/2/3+ {c['entrapment_1pep']:>3} {c['entrapment_2pep']:>2} {c['entrapment_3pep_plus']:>2}  "
                  f"({c['share_entrapment_1pep']:.0%} single)  T single {c['target_1pep']:>5} "
                  f"({c['share_target_1pep']:.1%})" + (f"  NO PEPTIDE {c['n_0pep']}" if c["n_0pep"] else ""))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {Path(a.out).name}: {len(rows)} rows")


if __name__ == "__main__":
    main()
