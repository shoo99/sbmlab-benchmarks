"""Known-answer test of make_entrapment.py, extract.py and score.py on synthetic data (run by CI; standard library).

1. Database: 200 random target proteins; a foreign proteome of 150 proteins, 10 of which carry a target peptide.
   foreign  -> the 10 are dropped, no entrapment peptide equals a target peptide, accessions carry ENTRAP_
   shuffle  -> every copy has its target's length, composition and K/R/P positions (so the same cleavage sites),
               shares no peptide with the targets, and r_peptides is about 1
2. Search outputs in the formats of DIA-NN (report.tsv), MaxQuant (combined/txt) and FragPipe (two experiment folders)
   with planted identifications: targets and entrapment hits with known q-values, decoys (MaxQuant 2.8's Decoy column
   in evidence and peptides, the older Reverse column in proteinGroups) and contaminants (not
   counted), a peptide shared by a target and an entrapment protein (counted as target), and a target-only search
   for the targets lost. extract.py's counts and score.py's estimates must equal the planted values exactly.
3. Per-run levels: DIA-NN rows per run at their run-level q; MaxQuant evidence with a second run of MBR transfers
   (MaxQuant 2.8's empty Type, and MULTI-MATCH) plus a MULTI-MATCH-NOVERIFY row that is left out; FragPipe
   combined_ion.tsv with MS/MS, MBR, MBC and unmatched ions. Transfers must not change the experiment-wide counts.
Usage: python3 test/test_synthetic.py
"""
import csv
import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
from make_entrapment import peptides, read_fasta, tryptic  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
rng = random.Random(11)
checks = []


def check(name, ok):
    checks.append((name, bool(ok)))


def protein(n):
    return "M" + "".join(rng.choice(AA) for _ in range(n - 1))


def write_tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


with tempfile.TemporaryDirectory() as tmp:
    t = Path(tmp)
    targets = [(f"sp|T{i:05d}|TGT{i}_HUMAN Target protein {i}", protein(rng.randint(200, 600))) for i in range(200)]
    with open(t / "target.fasta", "w") as fh:
        fh.writelines(f">{h}\n{s}\n" for h, s in targets)
    target_peps = set().union(*(peptides(s) for _, s in targets))
    foreign = []
    for i in range(150):
        s = protein(rng.randint(200, 600))
        if i < 10:                                            # plant a target peptide
            p = sorted(peptides(targets[i][1]))[0]
            s = s[:50] + "K" + p + s[50:]
        foreign.append((f"sp|F{i:05d}|FOR{i}_ARATH Foreign protein {i}", s))
    with open(t / "foreign.fasta", "w") as fh:
        fh.writelines(f">{h}\n{s}\n" for h, s in foreign)
    planted_drop = sum(1 for _, s in foreign if peptides(s) & target_peps)

    # ---- 1. databases
    run = lambda *a: subprocess.run([sys.executable, *a], check=True, capture_output=True, text=True)
    run(str(HERE / "make_entrapment.py"), str(t / "target.fasta"), "foreign", str(t / "foreign.fasta"), "--out", str(t / "db_foreign"))
    info = json.loads((t / "db_foreign" / "entrapment_info.json").read_text())
    ent = read_fasta(t / "db_foreign" / "entrapment.fasta")
    check(f"foreign: {planted_drop} proteins sharing a target peptide dropped", info["dropped_proteins"] == planted_drop
          and info["entrapment_proteins"] == 150 - planted_drop == len(ent))
    check("foreign: no entrapment peptide equals a target peptide", info["shared_peptides"] == 0)
    check("foreign: every accession carries ENTRAP_", all(h.split("|")[1].startswith("ENTRAP_F") for h, _ in ent))
    check("foreign: combined database = targets + entrapment",
          len(read_fasta(t / "db_foreign" / "target_plus_entrapment.fasta")) == 200 + len(ent))

    run(str(HERE / "make_entrapment.py"), str(t / "target.fasta"), "shuffle", "1", "--out", str(t / "db"))
    info = json.loads((t / "db" / "entrapment_info.json").read_text())
    ent = dict((h.split("|")[1].replace("ENTRAP_", ""), s) for h, s in read_fasta(t / "db" / "entrapment.fasta"))
    tgt = dict((h.split("|")[1], s) for h, s in targets)
    same = all(len(ent[a]) == len(tgt[a]) and sorted(ent[a]) == sorted(tgt[a])
               and [i for i, c in enumerate(ent[a]) if c in "KRP"] == [i for i, c in enumerate(tgt[a]) if c in "KRP"]
               for a in ent)
    check(f"shuffle: {len(ent)} copies keep length, composition and K/R/P positions", same and len(ent) >= 195)
    check("shuffle: same tryptic cleavage sites", all([len(p) for p in tryptic(ent[a])] == [len(p) for p in tryptic(tgt[a])] for a in ent))
    check("shuffle: no shared peptide with the targets", info["shared_peptides"] == 0)
    check(f"shuffle: r_peptides {info['r_peptides']} within 0.97-1.0", 0.97 <= info["r_peptides"] <= 1.0)

    # ---- 2. planted identifications
    tacc = [h.split(" ")[0] for h, _ in targets]
    eacc = [a.replace("|T", "|ENTRAP_T") for a in tacc]
    seqs = iter(sorted(target_peps))
    planted = []                          # (kind, peptide, proteins, q)
    for i in range(400):
        planted.append(("target", next(seqs), [tacc[i % 200]], round(rng.uniform(0, 0.05), 5)))
    for i in range(40):
        planted.append(("entrapment", next(seqs), [eacc[i]], round(rng.uniform(0, 0.05), 5)))
    planted.append(("target", next(seqs), [tacc[3], eacc[7]], 0.004))          # shared: counts as target
    for i in range(5):
        planted.append(("decoy", next(seqs), [f"rev_{tacc[i]}"], 0.003))
    for i in range(3):
        planted.append(("contaminant", next(seqs), [f"CON__P0{i}"], 0.002))
    grid = [0.001, 0.005, 0.01, 0.02, 0.05]
    want = {g: (sum(k == "target" and q <= g for k, _, _, q in planted), sum(k == "entrapment" and q <= g for k, _, _, q in planted))
            for g in grid}
    n_t_all = sum(k == "target" for k, *_ in planted); n_e_all = sum(k == "entrapment" for k, *_ in planted)

    def groups(key_of, thr=None):
        """Expected protein-level counts: one group per key, its lowest q, its kind (target wins over entrapment)."""
        g = {}
        for k, p, pr, q in planted:
            if k in ("decoy", "contaminant"):
                continue
            key = key_of(pr)
            q0, k0 = g.get(key, (1.0, "entrapment"))
            g[key] = (min(q0, q), "target" if "target" in (k0, k) else "entrapment")
        sel = [k for q, k in g.values() if thr is None or q <= thr]
        return sel.count("target"), sel.count("entrapment")

    # DIA-NN: two runs, Global.Q.Value is the planted q; protein groups one per peptide
    hdr = ["Run", "Precursor.Id", "Stripped.Sequence", "Protein.Ids", "Protein.Group", "Q.Value", "Global.Q.Value",
           "PG.Q.Value", "Global.PG.Q.Value"]
    rows = []
    for run_ in ("R1", "R2"):
        for k, pep, prots, q in planted:
            rows.append([run_, pep + "2", pep, ";".join(prots), ";".join(prots), min(0.05, q * 1.5), q, q, q])
    write_tsv(t / "diann" / "report.tsv", hdr, rows)
    base_rows = [r for r in rows if "ENTRAP_" not in r[3] or ";" in r[3]] + \
        [["R1", f"PEPTIDEXK{i}2", f"PEPTIDEXK{i}", tacc[i], tacc[i], 0.001, 0.001, 0.001, 0.001] for i in range(20)]
    write_tsv(t / "diann_base" / "report.tsv", hdr, base_rows)
    # MaxQuant: reported identifications only
    mq = t / "mq" / "combined" / "txt"
    flag = lambda k: ["+" if k == "decoy" else "", "+" if k == "contaminant" else ""]
    mbr_types = ["", "MULTI-MATCH"]                       # MaxQuant 2.8 leaves Type empty for MBR transfers
    write_tsv(mq / "evidence.txt", ["Raw file", "Modified sequence", "Charge", "Proteins", "Decoy", "Potential contaminant",
                                    "Type", "Match time difference"],
              [["R1", f"_{p}_", 2, ";".join(pr), *flag(k), "MULTI-MSMS", ""] for k, p, pr, q in planted]
              + [["R2", f"_{p}_", 2, ";".join(pr), *flag(k), mbr_types[i % 2], "0.01"] for i, (k, p, pr, q) in enumerate(planted[:100])]
              + [["R3", f"_{p}_", 2, ";".join(pr), *flag(k), "MULTI-MATCH-NOVERIFY", "0.01"] for k, p, pr, q in planted[:5]])
    mq_mbr = (sum(k == "target" for k, *_ in planted[:100]), sum(k == "entrapment" for k, *_ in planted[:100]))
    write_tsv(mq / "peptides.txt", ["Sequence", "Proteins", "Decoy", "Potential contaminant"],
              [[p, ";".join(pr), "+" if k == "decoy" else "", "+" if k == "contaminant" else ""] for k, p, pr, q in planted])
    write_tsv(mq / "proteinGroups.txt", ["Protein IDs", "Q-value", "Reverse", "Potential contaminant", "Only identified by site"],
              [[";".join(pr), q, "+" if k == "decoy" else "", "+" if k == "contaminant" else "", ""] for k, p, pr, q in planted])
    # FragPipe: two experiments, each with half of the identifications plus 50 in both
    for e, part in (("exp1", planted[:250]), ("exp2", planted[200:])):
        d = t / "fp" / e
        write_tsv(d / "psm.tsv", ["Peptide", "Modified Peptide", "Charge", "Protein", "Mapped Proteins"],
                  [[p, "", 2, pr[0], ", ".join(pr[1:])] for k, p, pr, q in part])
        write_tsv(d / "peptide.tsv", ["Peptide", "Protein", "Mapped Proteins"], [[p, pr[0], ", ".join(pr[1:])] for k, p, pr, q in part])
        write_tsv(d / "protein.tsv", ["Protein", "Indistinguishable Proteins"], [[pr[0], ", ".join(pr[1:])] for k, p, pr, q in part])
    match = ["MBR", "MBC", "unmatched"]
    write_tsv(t / "fp" / "combined_ion.tsv", ["Peptide Sequence", "Modified Sequence", "Charge", "Protein", "Mapped Proteins",
                                              "S1 Match Type", "S2 Match Type"],
              [[p, "", 2, pr[0], ", ".join(pr[1:]), "MS/MS", match[i % 3]] for i, (k, p, pr, q) in enumerate(planted)])
    kinds_by = lambda m: (sum(k == "target" for i, (k, *_) in enumerate(planted) if match[i % 3] in m),
                          sum(k == "entrapment" for i, (k, *_) in enumerate(planted) if match[i % 3] in m))
    write_tsv(t / "searches.tsv", ["search", "engine", "version", "mbr", "result", "baseline", "entrapment"],
              [["DIA-NN synthetic", "diann", "DIA-NN test", "on", "diann/report.tsv", "diann_base/report.tsv", "db/entrapment_info.json"],
               ["MaxQuant synthetic", "maxquant", "MaxQuant test", "off", "mq/combined/txt", "", "db/entrapment_info.json"],
               ["FragPipe synthetic", "fragpipe", "FragPipe test", "on", "fp", "", "db/entrapment_info.json"]])
    run(str(HERE / "extract.py"), str(t / "searches.tsv"), "--out", str(t / "counts.csv"))
    run(str(HERE / "score.py"), str(t / "counts.csv"), str(t / "out"))
    counts = list(csv.DictReader(open(t / "counts.csv")))
    get = lambda s, lev, th: next(c for c in counts if c["search"] == s and c["level"] == lev and c["threshold"] == th)
    for g in grid:
        c = get("DIA-NN synthetic", "precursor", str(g))
        check(f"DIA-NN precursors q <= {g}: target {want[g][0]}, entrapment {want[g][1]}",
              (int(c["n_target"]), int(c["n_entrapment"])) == want[g])
    c = get("DIA-NN synthetic", "protein", "0.01")
    check(f"DIA-NN protein groups at q <= 0.01: {groups(';'.join, 0.01)}",
          (int(c["n_target"]), int(c["n_entrapment"])) == groups(";".join, 0.01))
    c = get("DIA-NN synthetic", "precursor", "0.01")
    base_t = sum(1 for r in base_rows if r[0] == "R1" and float(r[6]) <= 0.01 and not r[3].startswith(("rev_", "CON__")))
    check(f"DIA-NN target-only search: {base_t} targets at q <= 0.01", int(c["n_baseline_target"]) == base_t)
    for s in ("MaxQuant synthetic", "FragPipe synthetic"):
        for lev in ("precursor", "peptide"):
            c = get(s, lev, "reported")
            check(f"{s} {lev}: target {n_t_all}, entrapment {n_e_all} (decoys and contaminants left out)",
                  (int(c["n_target"]), int(c["n_entrapment"])) == (n_t_all, n_e_all))
    c = get("MaxQuant synthetic", "protein", "0.01")
    check(f"MaxQuant protein groups use the Q-value column: {groups(';'.join, 0.01)} at q <= 0.01",
          (int(c["n_target"]), int(c["n_entrapment"])) == groups(";".join, 0.01))
    c = get("FragPipe synthetic", "protein", "reported")
    check(f"FragPipe proteins: {groups(lambda pr: pr[0])}", (int(c["n_target"]), int(c["n_entrapment"])) == groups(lambda pr: pr[0]))

    for g in grid:
        rq = lambda q: min(0.05, q * 1.5) <= g
        w = (2 * sum(k == "target" and rq(q) for k, _, _, q in planted), 2 * sum(k == "entrapment" and rq(q) for k, _, _, q in planted))
        c = get("DIA-NN synthetic", "precursor_run", str(g))
        check(f"DIA-NN per-run precursors at run-level q <= {g}: {w}", (int(c["n_target"]), int(c["n_entrapment"])) == w)
    check("DIA-NN: no transfer level (DIA-NN does not mark transfers)",
          not any(c["search"] == "DIA-NN synthetic" and c["level"] == "precursor_run_mbr" for c in counts))
    c = get("MaxQuant synthetic", "precursor_run", "reported")
    w = (n_t_all + mq_mbr[0], n_e_all + mq_mbr[1])
    check(f"MaxQuant per-run precursors: MS/MS + MBR = {w}, NOVERIFY left out", (int(c["n_target"]), int(c["n_entrapment"])) == w)
    c = get("MaxQuant synthetic", "precursor_run_mbr", "reported")
    check(f"MaxQuant transfers (empty Type and MULTI-MATCH): {mq_mbr}", (int(c["n_target"]), int(c["n_entrapment"])) == mq_mbr)
    c = get("FragPipe synthetic", "precursor_run", "reported")
    w = tuple(a + b for a, b in zip((n_t_all, n_e_all), kinds_by(("MBR", "MBC"))))
    check(f"FragPipe per-run ions: MS/MS + MBR + MBC = {w}, unmatched left out", (int(c["n_target"]), int(c["n_entrapment"])) == w)
    c = get("FragPipe synthetic", "precursor_run_mbr", "reported")
    check(f"FragPipe transfers (MBR): {kinds_by(('MBR',))}", (int(c["n_target"]), int(c["n_entrapment"])) == kinds_by(("MBR",)))

    scores = list(csv.DictReader(open(t / "out" / "scores.csv")))
    sv = lambda tool, sub, setting, metric: float(next(r["value"] for r in scores if r["tool"] == tool and r["subset"] == sub
                                                       and r["setting"] == setting and r["metric"] == metric))
    r_pep = info["r_peptides"]
    nt, ne = want[0.01]
    check("fdp_lower = N_E / (N_T + N_E)", abs(sv("DIA-NN synthetic", "precursor", "MBR on; q <= 0.01", "fdp_lower") - round(ne / (nt + ne), 5)) < 1e-9)
    check("fdp_combined = N_E (1 + 1/r) / (N_T + N_E)",
          abs(sv("DIA-NN synthetic", "precursor", "MBR on; q <= 0.01", "fdp_combined") - round(ne * (1 + 1 / r_pep) / (nt + ne), 5)) < 1e-9)
    check("share of targets lost = 1 - N_T / N_T(target-only)",
          abs(sv("DIA-NN synthetic", "precursor", "MBR on; q <= 0.01", "share_targets_lost") - round(1 - nt / base_t, 4)) < 1e-9)
    with open(t / "out" / "scores.csv") as fh:
        check("scores.csv has the standard header",
              next(csv.reader(fh)) == ["folder", "dataset", "tool", "version", "setting", "subset", "metric", "value", "n"])

for name, ok in checks:
    print(("ok   " if ok else "FAIL ") + name)
if not all(ok for _, ok in checks):
    sys.exit(1)
print(f"all {len(checks)} synthetic checks passed")
