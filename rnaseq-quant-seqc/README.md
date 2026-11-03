# RNA-seq quantifiers scored on SEQC/MAQC mixtures and ERCC spike-ins (in progress)

Tracking issue: [#4](https://github.com/shoo99/sbmlab-benchmarks/issues/4). No post yet.

Which RNA-seq quantifier gets known ratios right: Salmon, kallisto, STAR + featureCounts, or RSEM? This is the RNA-seq counterpart of the HYE proteomics benchmark.

**Status.** The pipeline and the scoring are finished and tested end to end on synthetic data with a known answer. That test runs in CI (`test/test_pipeline.sh`). The SEQC run itself has not been done yet, and the Results section will be filled in from it.

## Data

- SEQC/MAQC-III RNA-seq. The reads are on GEO/SRA: series GSE47792, SRA study SRP025982. Choose the runs with `./list_runs.sh`. The plan is one sequencing site and four replicates of each of A, B, C and D. The run accessions used go into `runs.tsv`.
- References:
  - GENCODE human release 46 (set `GENCODE` to change it): GRCh38 primary assembly, comprehensive CHR annotation, transcript sequences
  - ERCC92 sequences and the Mix 1 / Mix 2 concentration table from Thermo Fisher
- None of these files is in the repository; `download.sh` fetches them.

## Truth

- **Samples**:
  - A = Universal Human Reference RNA + ERCC Mix 1
  - B = Human Brain Reference RNA + ERCC Mix 2
  - C = 3 parts A + 1 part B, by total RNA
  - D = 1 part A + 3 parts B
- **ERCC**: the 92 spikes fall into four subgroups of 23, with Mix 1 / Mix 2 = 4, 1, 1/1.5 and 1/2. So the true log2(A/B) of each spike is 2, 0, −0.585 or −1.
  - Expression is per million *endogenous* reads, so the ERCC ratios shift by the difference in mRNA content between A and B.
  - Scores are therefore reported uncentred and centred on the median of the 1:1 subgroup. This is the analogue of centring HYE on the human median.
- **Titration**: a gene's expression in C and D must lie between A and B in the order A > C > D > B, or the reverse. This needs no knowledge of the true expression.
- **Mixture model**: because mixing was by total RNA, the mRNA fraction coming from A is 0.75w / (0.75w + 0.25) in C and 0.25w / (0.25w + 0.75) in D, where w is the mRNA content of A relative to B. `score_seqc.R` fits w and reports how far C and D are from the prediction.

## Tools

| Tool | How it is run (`run.sh`) |
|---|---|
| Salmon | index on transcripts + ERCC (`SALMON_DECOYS=1` for a decoy-aware index); `quant -l A` |
| kallisto | index on the same transcripts; `quant` (strandedness from `STRAND`) |
| STAR + featureCounts | STAR on genome + ERCC with the GTF; featureCounts per gene (`-p --countReadPairs`, `-s` from `STRAND`) |
| RSEM | `rsem-calculate-expression --alignments` on STAR's `Aligned.toTranscriptome.out.bam` |

Versions are recorded in `versions.tsv` and run times in `timings.tsv`. CI uses the Ubuntu 24.04 packages: Salmon 1.10.2, kallisto 0.48.0, STAR 2.7.11b, subread 2.0.6, and RSEM reporting 1.3.1. For the SEQC run, use current releases and record them.

## How to run

```bash
./list_runs.sh > all_runs.tsv        # choose runs, then list them in runs.tsv (sample, replicate, run_accession)
./download.sh ~/seqc 0               # 0 = all reads; or e.g. 10000000 read pairs per replicate (needs seqtk)
STRAND=none ./run.sh ~/seqc 24       # check work/salmon/*/lib_format_counts.json and set STRAND if the library is stranded
Rscript gene_matrix.R ~/seqc/work .
Rscript score_seqc.R . ~/seqc/ref/ercc_truth.tsv . "SEQC <site>" 3 > score_seqc.out
```

| Script | What it does | Output |
|---|---|---|
| `list_runs.sh` | Lists the study's runs at ENA | stdout |
| `download.sh` | Downloads reads (lanes concatenated, optional subsampling with seed 11), GENCODE and ERCC files | `<data>/fastq`, `<data>/ref` |
| `ercc_truth.py` | Converts Thermo's concentration table to `ercc_truth.tsv`, checking 92 spikes and the 4 / 1 / 0.67 / 0.5 design | `ref/ercc_truth.tsv` |
| `run.sh` | Runs the four quantifiers | `<data>/work`, with `versions.tsv` and `timings.tsv` |
| `gene_matrix.R` | Builds gene × sample counts per tool. Transcript estimates are summed per gene; ERCC spikes are their own genes | `genes_<tool>.tsv.gz` |
| `score_seqc.R` | Scores ERCC, titration order and the mixture model with `common/known_ratio.R` | printout, `scores.csv` |
| `simulate/simulate.py` | Builds a synthetic SEQC-like data set with a known answer: small spliced genome, isoforms, 16 ERCC-like spikes, w = 1.3, replicate noise | `ref/`, `fastq/`, `truth.json` |
| `test/test_pipeline.sh` | Runs simulate → run → matrix → score and checks that the simulated answer comes back (`test/check_synthetic.py`) | pass / fail |

On the synthetic data, all four tools recover the answer:
- fitted w 1.313–1.320 vs the true 1.3
- ERCC within ±0.5 after centring: 1.00
- titration order consistent for 0.89–0.91 of the changed genes

Two broken inputs fail the checks as they should: a wrong w in the truth file, and ERCC ratios with the sign flipped.

## Results

Pending the SEQC run.

## Traps

- kallisto 0.48 does not create missing parent directories for `-o`; `run.sh` creates them.
- STAR needs `--genomeSAindexNbases` scaled down for small genomes (the synthetic test); `run.sh` uses STAR's formula, capped at 14.
- GENCODE transcript headers carry `|`-separated fields; `run.sh` cuts them to the transcript ID so that all tools use the same names.
- featureCounts 2.0.2 and later needs `--countReadPairs` together with `-p` to count fragments rather than reads.

## Not included

FASTQ files, the GENCODE and ERCC references, and the `work/` directory with indices, alignments and per-tool outputs. The gene matrices (`genes_<tool>.tsv.gz`), `scores.csv`, `score_seqc.out`, `versions.tsv` and `timings.tsv` from the SEQC run will be committed.
