"""Builds an entrapment database for the HYE searches (issue #5; Wen et al. 2025, Nat Methods).

Any identification from an entrapment protein is false, because those proteins cannot be in the sample. Two kinds:
  foreign <fasta>   proteins of a species that cannot be in the sample (default plan: Arabidopsis thaliana Swiss-Prot).
                    Proteins sharing any tryptic peptide with the targets are dropped (conserved histones, ubiquitin, ...).
  shuffle <k>       k shuffled copies of every target protein. Within each tryptic peptide the residues other than K, R
                    and P are shuffled, so masses, cleavage sites and lengths stay the same; a shuffled peptide that
                    equals a target peptide is reshuffled (up to 20 tries, else that protein copy is dropped).
Peptides are compared after I -> L (isobaric). "Tryptic peptide" here: cleavage after K or R, not before P, no missed
cleavage, length 7-30 (the search settings of diann-hye). The entrapment-to-target ratio r is reported for proteins
and for these peptides; Wen et al. use the peptide ratio in the combined estimator.

Output (next to this script unless --out is given):
  <out>/entrapment.fasta           the entrapment proteins; accessions prefixed ENTRAP_ (sp|ENTRAP_P12345|...)
  <out>/target_plus_entrapment.fasta   targets followed by the entrapment proteins: the search database
  <out>/entrapment_info.json       mode, counts, r (protein and peptide), seed
Usage: python3 make_entrapment.py <target.fasta> foreign <foreign.fasta> [--out DIR]
       python3 make_entrapment.py <target.fasta> shuffle <k> [--out DIR] [--seed N]
Standard library only.
"""
import argparse
import json
import random
import re
from pathlib import Path

PREFIX = "ENTRAP_"
MIN_LEN, MAX_LEN = 7, 30


def read_fasta(path):
    out, header, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    out.append((header, "".join(seq)))
                header, seq = line[1:], []
            elif line:
                seq.append(line.strip())
    if header is not None:
        out.append((header, "".join(seq)))
    return out


def write_fasta(path, records):
    with open(path, "w") as fh:
        for h, s in records:
            fh.write(f">{h}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


def tryptic(seq):
    """Fully tryptic pieces: cleave after K/R unless followed by P (every piece, any length)."""
    return [p for p in re.split(r"(?<=[KR])(?!P)", seq) if p]


def peptides(seq):
    return {p.replace("I", "L") for p in tryptic(seq) if MIN_LEN <= len(p) <= MAX_LEN}


def entrap_header(h):
    """sp|P12345|NAME_HUMAN desc -> sp|ENTRAP_P12345|NAME_HUMAN desc; other headers: ENTRAP_<first word> ..."""
    m = re.match(r"^(sp|tr)\|([^|]+)\|(.*)$", h)
    if m:
        return f"{m.group(1)}|{PREFIX}{m.group(2)}|{m.group(3)}"
    return PREFIX + h


def shuffle_peptide(p, rng):
    """Shuffles the residues of p among the positions that are not K, R or P. K, R and P stay in place, so every
    cleavage site (K/R not before P) is kept, inside the peptide and at its boundaries, and the mass is unchanged."""
    pos = [i for i, c in enumerate(p) if c not in "KRP"]
    res = [p[i] for i in pos]
    rng.shuffle(res)
    q = list(p)
    for i, c in zip(pos, res):
        q[i] = c
    return "".join(q)


def shuffled_protein(seq, target_peps, rng, tries=20):
    """Shuffles each tryptic piece; a piece of counted length that equals a target peptide is reshuffled up to `tries`
    times. Returns None if some piece cannot be made different (for example a peptide of one repeated residue)."""
    pieces = []
    for p in tryptic(seq):
        q = shuffle_peptide(p, rng)
        if MIN_LEN <= len(p) <= MAX_LEN:
            n = 1
            while q.replace("I", "L") in target_peps and n < tries:
                q = shuffle_peptide(p, rng)
                n += 1
            if q.replace("I", "L") in target_peps:
                return None
        pieces.append(q)
    return "".join(pieces)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("target")
    ap.add_argument("mode", choices=["foreign", "shuffle"])
    ap.add_argument("arg", help="foreign FASTA, or number of shuffled copies k")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent))
    ap.add_argument("--seed", type=int, default=5)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    targets = read_fasta(a.target)
    target_peps = set()
    for _, s in targets:
        target_peps |= peptides(s)
    entrap, dropped = [], 0
    if a.mode == "foreign":
        for h, s in read_fasta(a.arg):
            if peptides(s) & target_peps:
                dropped += 1
                continue
            entrap.append((entrap_header(h), s))
    else:
        k = int(a.arg)
        rng = random.Random(a.seed)
        for copy in range(1, k + 1):
            for h, seq in targets:
                new = shuffled_protein(seq, target_peps, rng)
                if new is None:
                    dropped += 1
                    continue
                hh = entrap_header(h)
                if k > 1:
                    hh = hh.replace(PREFIX, f"{PREFIX}{copy}_", 1)
                entrap.append((hh, new))
    entrap_peps = set()
    for _, s in entrap:
        entrap_peps |= peptides(s)
    shared = len(entrap_peps & target_peps)
    info = {"mode": a.mode, "source": a.arg if a.mode == "foreign" else f"{a.arg} shuffled copies", "seed": a.seed,
            "target_proteins": len(targets), "entrapment_proteins": len(entrap), "dropped_proteins": dropped,
            "target_peptides": len(target_peps), "entrapment_peptides": len(entrap_peps),
            "shared_peptides": shared,
            "r_proteins": round(len(entrap) / len(targets), 4),
            "r_peptides": round(len(entrap_peps) / len(target_peps), 4),
            "peptide_rule": f"trypsin (K/R, not before P), no missed cleavage, {MIN_LEN}-{MAX_LEN} aa, I = L"}
    write_fasta(out / "entrapment.fasta", entrap)
    write_fasta(out / "target_plus_entrapment.fasta", targets + entrap)
    (out / "entrapment_info.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
