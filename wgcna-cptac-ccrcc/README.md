# WGCNA on CPTAC clear-cell RCC: independent re-run of the key numbers

Post: [Does WGCNA Really Work? Running It End to End on Public Proteomics Data](https://sbmlab.com/blog/does-wgcna-really-work-public-proteomics-end-to-end) · [한국어](https://sbmlab.com/blog/does-wgcna-really-work-public-proteomics-end-to-end-ko)

These scripts were written separately from the post's original analysis, from the data and settings stated in the post only, to check whether its main numbers come back. They do.

## Run

```bash
./download.sh   # ~48 MB: LinkedOmics CPTAC-CCRCC tumour proteome + clinical, TJ-RCC protein matrix (Zenodo 8063124, CC-BY-4.0); checks md5
./run_all.sh    # ~4.5 min on 16 threads; installs WGCNA into ./Rlib if it is missing
```

R 4.3.3, WGCNA 1.74 (the post used 1.73), impute 1.76.0. Output goes to `logs/`; the logs committed here are from the run on 2026-10-23.

## Settings

- Tumours only: 110 CPTAC tumours with grade G1 7 · G2 51 · G3 41 · G4 11.
- Proteins shared with TJ-RCC, ≤10% missing in both cohorts, **exact zeros in the TJ-RCC matrix treated as missing** (the matrix has no NA; 23% of entries are 0). 3,392 proteins pass; the top 3,000 by CPTAC variance are used.
- `blockwiseModules`: signed network, signed TOM, Pearson with pairwise-complete observations, power 10, minModuleSize 30, deepSplit 2, mergeCutHeight 0.25, one block, `pamRespectsDendro = FALSE`, randomSeed 12345.
- Grade is used as a number (1–4). Endothelial score = mean of z-scored PECAM1, VWF, CDH5, ENG, EMCN, ESAM from the full proteome.
- Random-set null: 1,000 random sets of the same size from the 3,000 input proteins; score = first principal component; compared on |r| with grade.
- AUC: G3–4 (52) vs G1–2 (58), logistic regression, stratified 5-fold CV repeated 50 times. The top-50 proteins are re-selected inside every training fold; PCs are refit inside every fold; module eigengenes come from the full-data network (no labels involved).

## Results (variant A = the post's input set)

| Check | Post | This run |
|---|---|---|
| Modules at power 10 / grey | 15 | 15 / 13% |
| First power reaching R² 0.80 / 0.85 / 0.90 | 10 / 18 / 20 | 10 / 17 / 19 scanning every power; 10 / 18 / 20 on the tutorial grid (1–10, then 12–20 by 2) |
| Power 18: mean connectivity / grey | ~2 / 43% | 2.08 / 42.6% |
| Modules with grade p < 0.05 | 7 | 7 |
| Strongest module with grade | r −0.46 | r −0.463 (186 proteins; hubs RASIP1, TNS2, CDH13 … PLVAP) |
| That module vs endothelial score | r 0.90 | 0.907 (0.899 with the 4 marker proteins removed from the module) |
| Grade r after adjusting for the score | −0.21 | −0.206 (partial p ≈ 0.031) |
| Significant modules beyond the random-set 95th percentile | 2 of 7 | 2 of 7 |
| CV AUC: top-50 score / 15 PCs / 15 eigengenes | 0.73 / 0.72 / 0.70 | 0.744 / 0.723 / 0.709 (top 50 selected on the full data instead: 0.821) |

## How much the input set matters

| Variant | Input proteins | Modules | Grade p < 0.05 | Beat random sets |
|---|---|---|---|---|
| A (post) | shared, TJ-RCC zeros = missing | 15 | 7 | 2 |
| B | shared, TJ-RCC zeros kept as values | 12 | 3 | 3 |
| C | CPTAC only, ≤10% missing | 10 | 2 | 2 |

The count of "significant" modules depends on preprocessing. In every variant the vascular module is the strongest, and only 2–3 modules beat size-matched random sets.

## Not re-run here

GO enrichment, module preservation (Zsummary), the subsampling table, the per-setting ARI bars, the RNA comparison and the tumour + normal run (TJ-RCC, 100 + 50).
