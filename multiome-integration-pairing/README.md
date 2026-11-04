# Unpaired RNA + ATAC integration scored with 10x Multiome pairing (in progress)

Tracking issue: [#16](https://github.com/shoo99/sbmlab-benchmarks/issues/16). No post yet.

**Status: running.** The data are the 10x PBMC 10k Multiome, where RNA and ATAC come from the same nuclei, so every true pairing is known. The pairing is hidden from each method and used only for scoring. First results:

| Method | FOSCTTM | Label accuracy |
|---|---|---|
| naive baseline | 0.119 | 0.832 |
| Seurat CCA | 0.062 | 0.916 |

GLUE is training. LIGER needs a workstation, because the rliger package is not on conda and CRAN is not reachable from the cloud environment.

**Scripts**:

- `download.sh` fetches the inputs.
- `seurat_methods.R` runs the baseline and Seurat CCA.
- `run_glue.py` runs GLUE.
- `score.py` writes the scores, and `figure.py` draws the figure.
