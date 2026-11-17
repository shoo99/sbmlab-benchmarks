"""Paired entrapment estimate of the peptide- and protein-level FDP (Wen et al. 2025, Nat Methods, Eq. 4), for the shuffled entrapment.

  FDP* = (N_E + N_{E>=s>T} + 2 N_{E>T>=s}) / (N_T + N_E)

N_T, N_E: discovered target and entrapment peptides (as in extract.py: a peptide is entrapment when every protein it maps
to is ENTRAP_, decoys and contaminants are left out). Each entrapment peptide is paired with the target peptide at the
same position in the protein it was shuffled from (make_entrapment.py shuffle 1 keeps lengths and K/R/P positions, so the
positions line up). N_{E>=s>T}: discovered entrapment peptides whose paired target was not discovered. N_{E>T>=s}:
discovered entrapment peptides whose paired target was discovered with a worse score. The paired method needs r = 1
(one shuffled copy per target); it is an upper bound like the combined estimate, but tighter.
Scores (higher is better): DIA-NN -PEP (lowest PEP over the peptide's precursors), MaxQuant -PEP (peptides.txt),
FragPipe Probability (highest over the experiments' peptide.tsv; MBR-only rows have none). Ties (equal scores) are not counted as E > T in FDP*;
FDP*_ties counts them as E > T, an upper variant. Discovered = as in extract.py: DIA-NN peptide q <= 0.01 (lowest
Global.Q.Value of its precursors), MaxQuant and FragPipe as reported.
Usage: python3 paired.py <target.fasta> <entrapment.fasta> <searches.tsv> [--out paired.csv]
       (searches.tsv as for extract.py; only peptide level, only rows whose entrapment_info.json says mode shuffle)
"""
import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract import protein_class, rows_tsv, split_proteins  # noqa: E402
from make_entrapment import read_fasta  # noqa: E402

PREFIX = "ENTRAP_"


def acc(header):
    parts = header.split("|")
    return parts[1] if len(parts) > 2 else header.split()[0]


def peptides_diann(path, q=0.01):
    import pandas as pd
    d = pd.read_parquet(path, columns=["Stripped.Sequence", "Protein.Ids", "Global.Q.Value", "PEP", "Decoy"]) \
        if str(path).endswith(".parquet") else pd.read_csv(path, sep="\t")
    d = d[d["Decoy"].astype(str).isin(["0", "0.0", "False"])]
    g = d.groupby("Stripped.Sequence").agg(q=("Global.Q.Value", "min"), pep=("PEP", "min"), prot=("Protein.Ids", "first"))
    g = g[g.q <= q]
    return {s: (-float(r.pep), split_proteins(r.prot)) for s, r in g.iterrows()}


def peptides_maxquant(folder):
    out = {}
    for r in rows_tsv(Path(folder) / "peptides.txt"):
        if "+" in (r.get("Reverse", ""), r.get("Decoy", ""), r.get("Potential contaminant", "")):
            continue
        out[r["Sequence"]] = (-float(r["PEP"]) if r.get("PEP") not in (None, "", "NaN") else float("-inf"),
                              split_proteins(r["Proteins"]))
    return out


def peptides_fragpipe(folder):
    out = {}
    for f in sorted(Path(folder).rglob("peptide.tsv")):
        for r in rows_tsv(f):
            if not r.get("Probability"):        # MBR-only rows (no MS/MS in this experiment) carry no probability
                sc = float("-inf")
            else:
                sc = float(r["Probability"])
            prots = [r["Protein"]] + split_proteins(r.get("Mapped Proteins", ""))
            if r["Peptide"] not in out or sc > out[r["Peptide"]][0]:
                out[r["Peptide"]] = (sc, prots)
    return out


def proteins_diann(path, q=0.01):
    import pandas as pd
    d = pd.read_parquet(path, columns=["Protein.Group", "Global.PG.Q.Value", "PG.PEP", "Decoy"]) \
        if str(path).endswith(".parquet") else pd.read_csv(path, sep="\t")
    d = d[d["Decoy"].astype(str).isin(["0", "0.0", "False"])]
    g = d.groupby("Protein.Group").agg(q=("Global.PG.Q.Value", "min"), pep=("PG.PEP", "min"))
    g = g[g.q <= q]
    return {pg: (-float(r.pep), split_proteins(pg)) for pg, r in g.iterrows()}


def proteins_maxquant(folder, q=0.01):
    out = {}
    for r in rows_tsv(Path(folder) / "proteinGroups.txt"):
        if "+" in (r.get("Reverse", ""), r.get("Decoy", ""), r.get("Potential contaminant", "")) or r.get("Only identified by site") == "+":
            continue
        if r.get("Q-value") not in (None, "") and float(r["Q-value"]) > q:
            continue
        out[r["Protein IDs"]] = (float(r["Score"]), split_proteins(r["Protein IDs"]))
    return out


def proteins_fragpipe(folder):
    out = {}
    for f in sorted(Path(folder).rglob("protein.tsv")):
        for r in rows_tsv(f):
            sc = float(r["Protein Probability"]) if r.get("Protein Probability") else float("-inf")
            prots = [r["Protein"]] + split_proteins(r.get("Indistinguishable Proteins", ""))
            if r["Protein"] not in out or sc > out[r["Protein"]][0]:
                out[r["Protein"]] = (sc, prots)
    return out


READ = {"diann": peptides_diann, "maxquant": peptides_maxquant, "fragpipe": peptides_fragpipe}
READ_PROT = {"diann": proteins_diann, "maxquant": proteins_maxquant, "fragpipe": proteins_fragpipe}


def paired_counts_protein(groups):
    """groups: {key: (score, accessions)}. An entrapment group (all accessions ENTRAP_) is paired with the target
    accession it was shuffled from; the paired target counts as discovered when any discovered target group contains it."""
    disc_t_acc, disc_e = {}, {}
    for k, (sc, prots) in groups.items():
        c = protein_class(prots)
        accs = [acc(p) if "|" in p else p for p in prots]
        if c == "target":
            for a_ in accs:
                if not a_.startswith(PREFIX):
                    disc_t_acc[a_] = max(sc, disc_t_acc.get(a_, float("-inf")))
        elif c == "entrapment":
            disc_e[k] = (sc, accs)
    nt = sum(1 for k, (sc, prots) in groups.items() if protein_class(prots) == "target")
    ne = len(disc_e)
    absent = gt = ties = 0
    for k, (sc, accs) in disc_e.items():
        partners = [a_.replace(PREFIX, "", 1) for a_ in accs]
        found = [disc_t_acc[t] for t in partners if t in disc_t_acc]
        if not found:
            absent += 1
        elif sc > max(found):
            gt += 1
        elif sc == max(found):
            ties += 1
    n = nt + ne
    return {"n_target": nt, "n_entrapment": ne, "n_E_ge_s_gt_T": absent, "n_E_gt_T_ge_s": gt, "n_ties": ties, "n_unpaired": 0,
            "fdp_paired": round((ne + absent + 2 * gt) / n, 5) if n else 0.0,
            "fdp_paired_ties": round((ne + absent + 2 * (gt + ties)) / n, 5) if n else 0.0,
            "fdp_lower": round(ne / n, 5) if n else 0.0}


def paired_counts(peps, target_seq, entrap_seq):
    """peps: {peptide: (score, proteins)}. Returns the counts of Eq. 4 plus diagnostics."""
    disc_t, disc_e = {}, {}
    for p, (sc, prots) in peps.items():
        c = protein_class(prots)
        if c == "target":
            disc_t[p] = sc
        elif c == "entrapment":
            disc_e[p] = (sc, prots)
    n_unpaired_target_absent = n_e_gt_t = ties = unpaired = 0
    for p, (sc, prots) in disc_e.items():
        pair = None
        for pr in prots:
            a = acc(pr) if "|" in pr else pr
            src = a.replace(PREFIX, "", 1)
            es, ts = entrap_seq.get(a), target_seq.get(src)
            if es is None or ts is None:
                continue
            i = es.find(p)
            if i >= 0 and len(es) == len(ts):
                pair = ts[i:i + len(p)]
                break
        if pair is None:
            unpaired += 1
            n_unpaired_target_absent += 1          # conservative: treat as paired target not discovered
            continue
        if pair not in disc_t:
            n_unpaired_target_absent += 1
        elif sc > disc_t[pair]:
            n_e_gt_t += 1
        elif sc == disc_t[pair]:
            ties += 1
    nt, ne = len(disc_t), len(disc_e)
    n = nt + ne
    return {"n_target": nt, "n_entrapment": ne, "n_E_ge_s_gt_T": n_unpaired_target_absent, "n_E_gt_T_ge_s": n_e_gt_t,
            "n_ties": ties, "n_unpaired": unpaired,
            "fdp_paired": round((ne + n_unpaired_target_absent + 2 * n_e_gt_t) / n, 5) if n else 0.0,
            "fdp_paired_ties": round((ne + n_unpaired_target_absent + 2 * (n_e_gt_t + ties)) / n, 5) if n else 0.0,
            "fdp_lower": round(ne / n, 5) if n else 0.0}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("target"); ap.add_argument("entrapment"); ap.add_argument("manifest")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "paired.csv"))
    a = ap.parse_args()
    target_seq = {acc(h): s for h, s in read_fasta(a.target)}
    entrap_seq = {acc(h): s for h, s in read_fasta(a.entrapment)}
    base = Path(a.manifest).resolve().parent
    rows = []
    for s in rows_tsv(a.manifest):
        info = json.loads((base / s["entrapment"]).read_text())
        if info.get("mode") != "shuffle":
            continue
        for level, r, c in (("peptide", info["r_peptides"], paired_counts(READ[s["engine"]](base / s["result"]), target_seq, entrap_seq)),
                            ("protein", info["r_proteins"], paired_counts_protein(READ_PROT[s["engine"]](base / s["result"])))):
            n = c["n_target"] + c["n_entrapment"]
            c["fdp_combined"] = round(c["n_entrapment"] * (1 + 1 / r) / n, 5) if n else 0.0
            rows.append({"search": s["search"], "engine": s["engine"], "mbr": s["mbr"], "level": level, "r": r, **c})
            print(f"{s['search'][:44]:<44} {level:<7} T {c['n_target']:>6} E {c['n_entrapment']:>4} E>=s>T {c['n_E_ge_s_gt_T']:>4} "
                  f"E>T>=s {c['n_E_gt_T_ge_s']:>3} ties {c['n_ties']:>3} unp {c['n_unpaired']:>2}  lower {c['fdp_lower']:.4f} "
                  f"paired {c['fdp_paired']:.4f} (ties {c['fdp_paired_ties']:.4f}) combined {c['fdp_combined']:.4f}")
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"wrote {Path(a.out).name}: {len(rows)} rows")


if __name__ == "__main__":
    main()
