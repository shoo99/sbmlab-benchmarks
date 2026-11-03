"""Synthetic SEQC-like data with a known answer, for testing the pipeline end to end (standard library only).

Writes, under <out dir>:
  ref/genome.fa, ref/annotation.gtf, ref/transcripts.fa   small genome with spliced genes and GENCODE-style headers
  ref/ercc.fa, ref/ercc.gtf, ref/ercc_truth.tsv            16 ERCC-like spikes in subgroups A-D (Mix1/Mix2 = 4, 1, 1/1.5, 1/2)
  fastq/<A|B|C|D>_<rep>_R1.fastq.gz, _R2.fastq.gz          paired-end reads, unstranded, 0.2% substitution errors
  truth.json                                               w (mRNA content of A relative to B) and the mixture fractions

The design follows SEQC: A and B carry ERCC Mix 1 and Mix 2; C = 3 parts A + 1 part B and D = 1 + 3 by total RNA,
so the mRNA fraction coming from A is 0.75w / (0.75w + 0.25) in C and 0.25w / (0.25w + 0.75) in D.

Usage: python3 simulate.py <out dir> [--pairs 40000] [--seed 1] [--w 1.3]
"""
import argparse
import gzip
import json
import math
import random
from pathlib import Path

READ_LEN = 75
COMP = str.maketrans("ACGT", "TGCA")
SUBGROUP_LOG2 = {"A": 2.0, "B": 0.0, "C": math.log2(1 / 1.5), "D": -1.0}   # log2(Mix 1 / Mix 2)


def revcomp(s):
    return s.translate(COMP)[::-1]


def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def build_reference(rng, n_genes):
    chroms, gtf, tx = {}, [], []          # tx: (transcript_id, gene_id, sequence)
    genes_per_chrom = math.ceil(n_genes / 4)
    g = 0
    for c in range(1, 5):
        name, seq = f"chrS{c}", [rng.choice("ACGT") for _ in range(300)]
        for _ in range(genes_per_chrom):
            if g == n_genes:
                break
            g += 1
            gid = f"ENSGS{g:05d}.1"
            strand = rng.choice("+-")
            exons = []
            for e in range(rng.randint(3, 5)):
                if e:
                    seq += [rng.choice("ACGT") for _ in range(rng.randint(80, 400))]        # intron
                start = len(seq) + 1
                seq += [rng.choice("ACGT") for _ in range(rng.randint(100, 300))]
                exons.append((start, len(seq)))
            gtf.append(f'{name}\tSIM\tgene\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\tgene_id "{gid}"; gene_type "protein_coding"; gene_name "GENE{g}";')
            isoforms = [list(range(len(exons)))]
            if len(exons) > 3 and rng.random() < 0.6:
                skip = rng.randrange(1, len(exons) - 1)
                isoforms.append([i for i in range(len(exons)) if i != skip])
            for k, keep in enumerate(isoforms, 1):
                tid = f"ENSTS{g:05d}{k}.1"
                ex = [exons[i] for i in keep]
                gtf.append(f'{name}\tSIM\ttranscript\t{ex[0][0]}\t{ex[-1][1]}\t.\t{strand}\t.\tgene_id "{gid}"; transcript_id "{tid}"; gene_type "protein_coding";')
                for s, e in ex:
                    gtf.append(f'{name}\tSIM\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id "{gid}"; transcript_id "{tid}"; gene_type "protein_coding";')
                t = "".join("".join(seq[s - 1:e]) for s, e in ex)
                tx.append((tid, gid, t if strand == "+" else revcomp(t)))
            seq += [rng.choice("ACGT") for _ in range(rng.randint(200, 600))]          # intergenic
        chroms[name] = "".join(seq)
    return chroms, gtf, tx


def build_ercc(rng):
    spikes = []
    for k in range(16):
        sub = "ABCD"[k % 4]
        sid = f"ERCC-{k + 1:05d}"
        spikes.append((sid, f"DQ{k + 1:06d}", sub, "".join(rng.choice("ACGT") for _ in range(rng.randint(500, 1000)))))
    return spikes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--pairs", type=int, default=40000, help="read pairs per sample")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--w", type=float, default=1.3, help="mRNA content of A relative to B")
    ap.add_argument("--genes", type=int, default=120)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    out = Path(args.out)
    (out / "ref").mkdir(parents=True, exist_ok=True)
    (out / "fastq").mkdir(parents=True, exist_ok=True)

    chroms, gtf, tx = build_reference(rng, args.genes)
    spikes = build_ercc(rng)
    with open(out / "ref/genome.fa", "w") as f:
        f.writelines(f">{c}\n{wrap(s)}\n" for c, s in chroms.items())
    (out / "ref/annotation.gtf").write_text("\n".join(gtf) + "\n")
    with open(out / "ref/transcripts.fa", "w") as f:
        f.writelines(f">{t}|{g}|-|-|{t.split('.')[0]}-201|GENE|{len(s)}|protein_coding|\n{wrap(s)}\n" for t, g, s in tx)
    with open(out / "ref/ercc.fa", "w") as f:
        f.writelines(f">{sid}\n{wrap(s)}\n" for sid, _, _, s in spikes)
    (out / "ref/ercc.gtf").write_text("".join(
        f'{sid}\tERCC\texon\t1\t{len(s)}\t.\t+\t.\tgene_id "{sid}"; transcript_id "{dq}";\n' for sid, dq, _, s in spikes))
    with open(out / "ref/ercc_truth.tsv", "w") as f:
        f.write("ercc_id\tsubgroup\tlog2_mix1_over_mix2\n")
        f.writelines(f"{sid}\t{sub}\t{SUBGROUP_LOG2[sub]:.6f}\n" for sid, _, sub, _ in spikes)

    # mRNA composition of A and B (molecule fractions); half the genes change between A and B by 1-3 log2 units
    gene_of = {}
    for t, g, _ in tx:
        gene_of.setdefault(g, []).append(t)
    a, b = {}, {}
    for g, ts in gene_of.items():
        base = rng.lognormvariate(0, 1.2)
        d = rng.choice([-3, -2, -1, 1, 2, 3]) * rng.uniform(0.8, 1.2) if rng.random() < 0.5 else rng.gauss(0, 0.1)
        for t in ts:
            share = rng.uniform(0.3, 1.0)
            a[t], b[t] = base * share, base * share * 2 ** d
    sa, sb = sum(a.values()), sum(b.values())
    a = {t: v / sa for t, v in a.items()}
    b = {t: v / sb for t, v in b.items()}

    w, s = args.w, 0.04                                   # m_A / m_B, and ERCC amount relative to total RNA
    mix1 = {sid: 2 ** rng.uniform(-3, 3) for sid, _, _, _ in spikes}
    mix2 = {sid: mix1[sid] / 2 ** SUBGROUP_LOG2[sub] for sid, _, sub, _ in spikes}
    n1, n2 = sum(mix1.values()), sum(mix2.values())
    mix1 = {k: v / n1 for k, v in mix1.items()}
    mix2 = {k: v / n1 for k, v in mix2.items()}           # same scale as mix1: Mix 2 differs only by the ratios
    fA = {"A": 1.0, "B": 0.0, "C": 0.75, "D": 0.25}       # share of total RNA from A
    seqs = {t: q for t, _, q in tx} | {sid: q for sid, _, _, q in spikes}
    truth = {"w": w, "mrna_fraction_from_A": {}}
    for sample in "ABCD":
        f = fA[sample]
        mrna = f * w + (1 - f) * 1.0
        fa = f * w / mrna
        truth["mrna_fraction_from_A"][sample] = fa
        mol = {t: fa * a[t] + (1 - fa) * b[t] for t in a}                       # per unit mRNA
        for sid in mix1:
            mol[sid] = s * (f * mix1[sid] + (1 - f) * mix2[sid]) / mrna       # spikes per unit mRNA
        for rep in range(1, 5):
            names = list(mol)
            weights = [mol[t] * rng.lognormvariate(0, 0.15) * len(seqs[t]) for t in names]
            picks = rng.choices(names, weights=weights, k=args.pairs)
            with gzip.open(out / f"fastq/{sample}_{rep}_R1.fastq.gz", "wt", compresslevel=1) as r1, \
                 gzip.open(out / f"fastq/{sample}_{rep}_R2.fastq.gz", "wt", compresslevel=1) as r2:
                for i, t in enumerate(picks):
                    q = seqs[t]
                    flen = min(len(q), max(READ_LEN + 25, int(rng.gauss(220, 30))))
                    st = rng.randint(0, len(q) - flen)
                    frag = q[st:st + flen]
                    if rng.random() < 0.5:                                        # unstranded library
                        frag = revcomp(frag)
                    reads = [frag[:READ_LEN], revcomp(frag[-READ_LEN:])]
                    for fh, rd in zip((r1, r2), reads):
                        rd = "".join(rng.choice("ACGT".replace(c, "")) if rng.random() < 0.002 else c for c in rd)
                        fh.write(f"@{sample}{rep}_{i}/{1 if fh is r1 else 2}\n{rd}\n+\n{'I' * len(rd)}\n")
    (out / "truth.json").write_text(json.dumps(truth, indent=1) + "\n")
    print(f"wrote {len(tx)} transcripts of {len(gene_of)} genes, {len(spikes)} spikes, 16 samples x {args.pairs} pairs to {out}")


if __name__ == "__main__":
    main()
