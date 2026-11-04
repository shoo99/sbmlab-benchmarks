# scRNA-seq doublet detection scored on genotype- and species-labelled doublets (2026-11-04)

Tracking issue: [#7](https://github.com/shoo99/sbmlab-benchmarks/issues/7). No post yet.

How many real doublets do scDblFinder, scds, DoubletFinder and Scrublet find, and how many singlets do they wrongly flag? Two truths label doublets independently of expression:

- **Genotype.** Kang et al. 2018 pooled PBMC from 8 donors in each lane and demultiplexed them by SNPs with demuxlet. Droplets with two genotypes are doublets.
- **Species.** A 10x mixture of human (HEK293T) and mouse (NIH3T3) cells. Droplets with many reads from both genomes are doublets.

**Status: ready for a workstation run.** The pipeline runs end to end on synthetic data in the layout of both downloads, and CI tests it in the pinned conda environment. The downloads come from GEO and 10x Genomics, which this repository's cloud environment cannot reach, so there are no real-data numbers yet.

## Data

| Data set | Droplets | Truth | What the truth can see |
|---|---|---|---|
| `kang_ctrl`, `kang_stim` | about 14,600 and 14,400 (GEO [GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583), batch 2: control and IFN-β lanes, 8 donors each) | demuxlet `multiplets`: `doublet` = doublet, `singlet` = singlet, `ambs` = unknown | doublets of two donors only; a same-donor doublet is called a singlet. With 8 donors, about 7/8 of doublets are visible. `prepare.R` computes the share from the donor mix as 1 − Σ p². |
| `hgmm_12k` | about 12,000 (10x Genomics, Cell Ranger 2.1.0, hg19 + mm10) | minor-species share of UMIs ≥ 0.2: doublet; ≤ 0.1: singlet; in between: unknown (ambient RNA gives a few percent) | cross-species doublets only, about half of all doublets: 2p(1 − p) for human share p |

**Every droplet is scored, but recall and precision count only droplets with a known truth.** A tool that flags an invisible doublet (same donor, or same species) is counted as wrong. The precisions here are therefore lower bounds; `datasets.tsv` records the visible share for each data set.

**Expected doublet rate**: 0.8% per 1,000 droplets, the 10x rule of thumb. That gives about 12% for the Kang lanes and 10% for hgmm_12k.

## Tools

The pinned conda environment `conda-linux-64.lock` (368 packages with checksums; `environment.yml` lists the main ones) has:

| Tool | Version | Settings |
|---|---|---|
| scDblFinder (default rate) | 1.24.10 | `scDblFinder(sce)`, which uses its own rate (1% per 1,000 cells) |
| scDblFinder (expected rate) | 1.24.10 | `dbr` = the expected rate |
| scds cxds / bcds / hybrid | 1.26.1 | `cxds_bcds_hybrid()`; scores only, no calls |
| DoubletFinder | 2.0.6 (Seurat 5.5.1) | log-normalisation, 2,000 variable genes, PCs 1–10, pN 0.25, pK from `paramSweep` + `find.pK`, nExp = expected rate × droplets; score = pANN |
| Scrublet | 0.2.3 | `expected_doublet_rate` = the expected rate, `scrub_doublets()` defaults; calls from its own threshold |

R 4.5.3 and Python 3.12. Seeds: `set.seed(7)` before each R tool, `random_state = 0` for Scrublet.

## Scoring

`score.py` uses the standard library only.

| Metric | Meaning |
|---|---|
| `auprc` | average precision of the score over droplets with known truth (tied scores form one step) |
| `doublet_share` | the chance level of `auprc`: known doublets / known droplets |
| `recall_at_expected`, `precision_at_expected` | the top k droplets by score are called doublets, with k = expected rate × all droplets |
| `recall_own_call`, `precision_own_call` | the tool's own doublet calls (scDblFinder, DoubletFinder, Scrublet) |

The issue also asks about the effect on clustering (clusters made of doublets). This is not measured yet.

## How to run

Build the environment once:

```bash
micromamba create -p ./env --file conda-linux-64.lock
```

Download the data:

```bash
bash download.sh
```

Run the study:

```bash
DBL_RSCRIPT=./env/bin/Rscript DBL_PYTHON=./env/bin/python3 bash run.sh
```

- **What `run.sh` does**:
  1. `prepare.R` writes each data set with its truth to `work/` (not tracked).
  2. `run_r_tools.R` and `run_scrublet.py` score every droplet.
  3. `collect.py` bundles the per-droplet truth and scores into `cell_scores.tsv.gz` and `datasets.tsv`.
  4. `score.py` writes `scores.csv`.
- **Printouts**: `prepare.out`, `tools.out` (with seconds per tool) and `score.out`.
- **Options**:
  - `DF_SWEEP=0` skips DoubletFinder's pK sweep (pK 0.09)
  - `EXPECTED_RATE` overrides the expected rate
- **After the run**: commit `cell_scores.tsv.gz`, `datasets.tsv`, `scores.csv` and the `.out` files, add the results to this README, and change `ci` to `reproduce` in `catalog.csv` with a `score.py` re-run in `ci/reproduce.sh`. `score.py` needs only the two bundled tables.
- **Test**: `bash test/test_pipeline.sh` (same `DBL_*` variables) runs the whole pipeline on synthetic data. It takes about 4 min.
  - The synthetic data are two Kang-layout lanes (4 donors, 3 cell types, 10% doublets, demuxlet-style labels including `ambs` and same-donor "singlets") and an hgmm-layout mixture with ambient RNA.
  - It checks 31 known answers:
    - the truth `prepare.R` derives equals the generator's labels
    - all seven tools reach AUPRC of at least twice the doublet share
    - scDblFinder and Scrublet find at least half the visible doublets at the expected rate
  - CI runs this test in a job that builds the environment from the lock file.

## Results

Pending the workstation run.

## Traps

- **Neither truth sees all doublets.** Same-donor and same-species doublets look like singlets to the truth, so a tool that finds them is penalised. Compare tools within a data set, not against an absolute precision.
- **The expected rate decides a lot.** Recall and precision at the expected rate depend on k, and so do DoubletFinder's nExp and scDblFinder's `dbr`. Kang et al. loaded more cells than usual, because demuxlet can resolve doublets between donors, so the 10x rule of thumb may be off. Compare the demuxlet doublet share in `prepare.out` with the expected rate.
- **Cross-species doublets are easy.** A human + mouse doublet has expression no single cell has. Genotype doublets between two cells of the same type in different donors are invisible to expression-based tools by construction. Expect the species truth to give higher scores.

## Files

| File | Content |
|---|---|
| `download.sh` | fetches GSE96583 batch 2 and 10x hgmm_12k |
| `prepare.R` | builds each data set with its truth (`work/<dataset>/`) |
| `run_r_tools.R` | scDblFinder (2 settings), scds (3 scores), DoubletFinder |
| `run_scrublet.py` | Scrublet |
| `collect.py` | bundles per-droplet truth and scores (`cell_scores.tsv.gz`, `datasets.tsv`) |
| `score.py` | AUPRC, recall and precision → `scores.csv` |
| `run.sh` | the whole study |
| `conda-linux-64.lock`, `environment.yml` | the pinned environment |
| `test/make_synthetic.R`, `test/check_synthetic.py`, `test/test_pipeline.sh` | end-to-end known-answer test (CI) |
