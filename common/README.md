# common: shared scoring code

One implementation of the scoring rules in [CONTRIBUTING.md](../CONTRIBUTING.md), in R (`known_ratio.R`, base R only) and Python (`known_ratio.py`, numpy only). New folders use it instead of copying the functions; existing folders keep their own code. An existing folder is moved over only when a re-run gives byte-identical outputs, so no published number changes silently.

```r
source(file.path(repo, "common", "known_ratio.R"))
```
```python
sys.path.insert(0, str(repo / "common")); import known_ratio as kr
```

| Function | What it does |
|---|---|
| `species_from_names(x, suffixes)` | Group label from UniProt entry names (`_HUMAN`, `_YEAST`, `_ECOLI` by default). A group naming none or several groups gets NA |
| `species_from_taxa(x, taxa)` | The same from `OX=` taxon IDs in FASTA headers (MaxQuant). IDs not in the map are ignored |
| `log2_ratio(m, a, b, min_valid = 2)` | Turns zeros, negatives and non-finite values into missing, then takes log2. Returns `ok` (at least `min_valid` values in each group) and `fc` = mean log2(b) − mean log2(a) |
| `group_log2_means(m, groups, min_valid = 2)` | Mean log2 per group of columns, NA below `min_valid` |
| `score_groups(fc, ok, group, truth, tol = 0.5, centre = NULL)` | For each group with a known truth: n, median, IQR (R quantile type 7) and the share within `tol` of the truth. With `centre` (e.g. `"Human"`) it also gives the share after subtracting that group's median. Values are rounded to 3 decimals |
| `call_outcome(called, fc, group, truth)`, `call_summary(outcome)` | Labels each call as right / wrong direction / false (truth 0) / not called / no truth. FDP = (false + wrong direction) / calls with a truth |
| `scores_rows(...)` (R), `write_scores(rows, path)` | Rows in the standard `scores.csv` format |

Constants: `HYE_TRUTH` (log2 B/A = 0 / −1 / +2), `HYE_NAME_SUFFIX`, `HYE_TAXON`, `SCORES_COLUMNS`.

## Tests

`Rscript common/tests/test_known_ratio.R` and `python3 common/tests/test_known_ratio.py` run in CI. Each checks two things:

- **A hand-worked toy matrix.** It covers zeros, 2-of-3 filtering, a group naming two species, a lookalike suffix (`_HUMANX`), centring, and a wrong-direction call.
- **Re-creation of the bundled DIA-NN tables.** The functions alone must give exactly `diann-hye/diann_scores.csv` and `diann-2.7-hye/diann_scores_270.csv`. This shows that both implementations agree with each other and with the rules used so far.
