"""Known-answer test of paired.py (paired entrapment estimate, Wen et al. 2025 Eq. 4) on synthetic data (run by CI).

400 random target proteins, one shuffled copy each (make_entrapment.py shuffle 1). Planted discoveries, the same in
the DIA-NN, MaxQuant and FragPipe formats, at peptide level (one peptide per protein) and protein level:
  300 target peptides/proteins; 20 entrapment ones: 10 whose paired target was not discovered, 6 that score better
  than their discovered paired target, 3 that score worse, 1 tie.
Expected: N_T = 300, N_E = 20, N_{E>=s>T} = 10, N_{E>T>=s} = 6, ties = 1,
FDP* = (20 + 10 + 2*6) / 320, FDP*_ties = (20 + 10 + 2*7) / 320, at both levels and for all three engines.
Usage: python3 test/test_paired.py   (needs pandas, as the DIA-NN reader does)
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
from make_entrapment import peptides, read_fasta  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
rng = random.Random(7)
checks = []


def check(name, ok):
    checks.append((name, bool(ok)))


def write_tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


with tempfile.TemporaryDirectory() as tmp:
    t = Path(tmp)
    targets = [(f"sp|T{i:05d}|TGT{i}_HUMAN Target {i}", "M" + "".join(rng.choice(AA) for _ in range(rng.randint(200, 400))))
               for i in range(400)]
    with open(t / "target.fasta", "w") as fh:
        fh.writelines(f">{h}\n{s}\n" for h, s in targets)
    subprocess.run([sys.executable, str(HERE / "make_entrapment.py"), str(t / "target.fasta"), "shuffle", "1", "--out", str(t / "db")],
                   check=True, capture_output=True, text=True)
    tseq = {h.split("|")[1]: s for h, s in targets}
    eseq = {h.split("|")[1]: s for h, s in read_fasta(t / "db" / "entrapment.fasta")}
    # one target peptide per protein, and its paired entrapment peptide at the same position
    pairs = []
    for a, s in tseq.items():
        e = eseq.get("ENTRAP_" + a)
        if e is None:
            continue
        for p in sorted(peptides(s)):
            i = s.find(p)
            if i >= 0 and s.count(p) == 1 and e.count(e[i:i + len(p)]) == 1:
                pairs.append((a, p, e[i:i + len(p)]))
                break
    check(f"enough paired proteins ({len(pairs)})", len(pairs) >= 330)
    rng.shuffle(pairs)
    # planted: 300 discovered targets; entrapment cases among them (paired target discovered) and outside (absent)
    disc_t = pairs[:300]
    absent = pairs[300:310]                 # entrapment discovered, paired target not discovered
    gt, lt, tie = pairs[:6], pairs[6:9], pairs[9:10]
    tscore = {a: rng.uniform(0.5, 0.9) for a, _, _ in disc_t}
    escore = {}
    for a, _, _ in gt: escore[a] = tscore[a] + 0.05
    for a, _, _ in lt: escore[a] = tscore[a] - 0.05
    for a, _, _ in tie: escore[a] = tscore[a]
    for a, _, _ in absent: escore[a] = rng.uniform(0.5, 0.9)
    ent = gt + lt + tie + absent
    acc_t = lambda a: f"sp|{a}|TGT"
    acc_e = lambda a: f"sp|ENTRAP_{a}|TGT"
    # DIA-NN: score = -PEP, so PEP = 1 - score; one row per peptide; protein group = its protein
    rows = [[p + "2", p, acc_t(a), acc_t(a), 0.001, 1 - tscore[a], 0.001, 1 - tscore[a], 0] for a, p, _ in disc_t] + \
           [[ep + "2", ep, acc_e(a), acc_e(a), 0.001, 1 - escore[a], 0.001, 1 - escore[a], 0] for a, _, ep in ent]
    write_tsv(t / "diann" / "report.tsv", ["Precursor.Id", "Stripped.Sequence", "Protein.Ids", "Protein.Group", "Global.Q.Value",
                                           "PEP", "Global.PG.Q.Value", "PG.PEP", "Decoy"], rows)
    # MaxQuant: peptides.txt (score = -PEP) and proteinGroups.txt (Score)
    mq = t / "mq" / "combined" / "txt"
    write_tsv(mq / "peptides.txt", ["Sequence", "Proteins", "PEP", "Decoy", "Potential contaminant"],
              [[p, acc_t(a), 1 - tscore[a], "", ""] for a, p, _ in disc_t] + [[ep, acc_e(a), 1 - escore[a], "", ""] for a, _, ep in ent])
    write_tsv(mq / "proteinGroups.txt", ["Protein IDs", "Q-value", "Score", "Decoy", "Potential contaminant", "Only identified by site"],
              [[acc_t(a), 0.001, tscore[a], "", "", ""] for a, _, _ in disc_t] + [[acc_e(a), 0.001, escore[a], "", "", ""] for a, _, _ in ent])
    # FragPipe: two experiments; the second repeats half the rows with lower scores and one MBR-only row without probability
    for exp, rows_t, rows_e, delta in (("e1", disc_t, ent, 0.0), ("e2", disc_t[:150], ent[:10], -0.2)):
        d = t / "fp" / exp
        write_tsv(d / "peptide.tsv", ["Peptide", "Probability", "Protein", "Mapped Proteins"],
                  [[p, tscore[a] + delta, acc_t(a), ""] for a, p, _ in rows_t] + [[ep, escore[a] + delta, acc_e(a), ""] for a, _, ep in rows_e]
                  + ([[disc_t[0][1], "", acc_t(disc_t[0][0]), ""]] if exp == "e2" else []))
        write_tsv(d / "protein.tsv", ["Protein", "Protein Probability", "Indistinguishable Proteins"],
                  [[acc_t(a), tscore[a] + delta, ""] for a, _, _ in rows_t] + [[acc_e(a), escore[a] + delta, ""] for a, _, _ in rows_e])
    write_tsv(t / "searches.tsv", ["search", "engine", "version", "mbr", "result", "baseline", "entrapment"],
              [["DIA-NN synthetic", "diann", "t", "off", "diann/report.tsv", "", "db/entrapment_info.json"],
               ["MaxQuant synthetic", "maxquant", "t", "off", "mq/combined/txt", "", "db/entrapment_info.json"],
               ["FragPipe synthetic", "fragpipe", "t", "off", "fp", "", "db/entrapment_info.json"]])
    info = json.loads((t / "db" / "entrapment_info.json").read_text())
    check("shuffled entrapment database (mode shuffle)", info["mode"] == "shuffle")
    subprocess.run([sys.executable, str(HERE / "paired.py"), str(t / "target.fasta"), str(t / "db" / "entrapment.fasta"),
                    str(t / "searches.tsv"), "--out", str(t / "paired.csv")], check=True, capture_output=True, text=True)
    got = list(csv.DictReader(open(t / "paired.csv")))
    want = {"n_target": 300, "n_entrapment": 20, "n_E_ge_s_gt_T": 10, "n_E_gt_T_ge_s": 6, "n_ties": 1}
    for r in got:
        check(f"{r['search']} {r['level']}: counts {want}", all(int(r[k]) == v for k, v in want.items()))
        check(f"{r['search']} {r['level']}: FDP* = 42/320", abs(float(r["fdp_paired"]) - round(42 / 320, 5)) < 1e-9)
        check(f"{r['search']} {r['level']}: FDP*_ties = 44/320", abs(float(r["fdp_paired_ties"]) - round(44 / 320, 5)) < 1e-9)
    check("6 result rows (3 engines x peptide, protein)", len(got) == 6)

for name, ok in checks:
    print(("ok   " if ok else "FAIL ") + name)
if not all(ok for _, ok in checks):
    sys.exit(1)
print(f"all {len(checks)} paired checks passed")
