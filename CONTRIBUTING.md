# Adding to this repository

Each folder answers one measured question, usually for one post. The rules below keep the folders comparable and re-runnable. `ci/lint.py` checks the mechanical ones on every pull request.

## What a folder is for

- **Score against a right answer where one exists**: known mixing ratios, spike-ins, a curated truth set (for example GIAB), a simulation with planted effects, or shuffled labels whose expected result is chance. When there is no right answer (agreement between tools, speed), say so in the README's Truth section.
- **Use public data only**, so anyone can re-run it. Avoid data behind an access committee.
- **Record what you ran**: tool versions, non-default settings, the machine and how long it took.
- **Write down the traps**: what failed or surprised you, and how you found out.

## Adding a folder

1. Pick a name of the form `<topic>-<dataset>`, lowercase with hyphens (for example `rnaseq-quant-seqc`). Don't rename a folder once a post, the preprint or a Zenodo release links to it; GitHub does not redirect moved folders.
2. Copy `docs/template/README.md` and `docs/template/run.sh` into it and fill them in.
3. Add a row to `catalog.csv` (columns below).
4. If any script in the folder can re-run from files in the repository alone, add it to `ci/reproduce.sh` and set `ci` to `reproduce`.
5. Run `python3 ci/lint.py` (and `bash ci/reproduce.sh` if you added re-runs) before opening the pull request.

## What goes in, and what stays out

- **In**: the scripts as run, their printed output (`<script>.out`), result tables, small inputs you may redistribute, and the README.
- **Out**:
  - raw data (link the source instead)
  - files whose terms forbid redistribution, such as KEGG API responses
  - single files over about 10 MB; put large intermediates in a Zenodo record and link it from the README

## Scripts

- Take input paths as arguments (or one variable at the top) and **write outputs next to the script**. Never write to an absolute path; `ci/lint.py` flags write calls with one. The only exception is `string-hubs/figure.py`, which writes the post's figure into the blog repository; it is listed in `ci/lint.py` with its reason.
- End each printout with the versions used: R `sessionInfo()` or `packageVersion()`, Python `module.__version__`.
- Set random seeds and say so in the README.
- Shell scripts:
  - Use `cd <dir> || exit 1`.
  - Capture an exit status in the same command: `rc=0; tool ... || rc=$?`. `$(date)` resets `$?`, so `echo "$(date) exit $?"` always prints 0.

## Scoring rules for the HYE benchmark (PXD028735)

Folders that score on the human / yeast / *E. coli* mixture use the same rules unless their README says otherwise:

- zero intensities become missing; values are log2-transformed
- a protein needs at least 2 of 3 values in each condition
- log2(B/A) is the mean of log2 B minus the mean of log2 A
- "within ±0.5" means within 0.5 of the truth (human 0, yeast −1, *E. coli* +2), uncentred unless stated
- species must be unambiguous: every protein in the group comes from one species
- a differential call is correct only in the right direction (yeast down, *E. coli* up); human calls and wrong-direction calls both count as false

If a folder centres on the human median, or uses another ratio or filter, the README says so next to the numbers.

## Results in a standard table (new folders)

New folders also write their headline numbers to `scores.csv` in long format, one row per number:

| Column | Meaning | Example |
|---|---|---|
| `folder` | this folder | `rnaseq-quant-seqc` |
| `dataset` | data and subset scored | `SEQC A vs B` |
| `tool` | tool or method | `salmon` |
| `version` | its version | `1.10.3` |
| `setting` | the setting that differs between rows | `default` |
| `subset` | group the number is for | `ERCC` |
| `metric` | what is measured | `within_0.5` |
| `value` | the number | `0.912` |
| `n` | how many items it is computed on | `92` |

`ci/lint.py` checks the header of every `scores.csv`. With one format across folders, a summary figure or table can be generated for the whole repository.

## catalog.csv

One row per analysis folder (not `docs` or `ci`):

| Column | Content |
|---|---|
| `folder` | folder name |
| `area` | for example `proteomics`, `transcriptomics`, `single-cell`, `statistics`, `networks-enrichment`, `databases`, `multi-omics` |
| `data` | data set and source |
| `truth` | what the result is scored against, or `none (...)` with what is measured instead |
| `tools` | main tools with versions, separated by `;` |
| `post` | the post (or preprint) link |
| `ci` | `reproduce` (re-run and compared in CI), `smoke` (run in CI, output not compared) or `no` (needs data or services outside the repository) |

## Checks

`.github/workflows/checks.yml` runs on every pull request and every push to `main`:

- `python3 ci/lint.py`:
  - `catalog.csv` matches the folders, and its `ci` column matches `ci/reproduce.sh`
  - every `.json` file is strict JSON
  - every `.py` and `.sh` file parses
  - no script writes to an absolute path
  - every `scores.csv` has the standard header
- `bash ci/reproduce.sh r` runs on Ubuntu 24.04 packages: R 4.3.3 and limma 3.58.1, the versions the bundled outputs were made with. It checks that every `.R` file parses, then re-runs the R scripts that need only repository files and compares their output with the bundled files.
- `bash ci/reproduce.sh python` does the same for the Python scripts.

Re-runs happen in a temporary copy, so they never touch the working tree. Printouts are compared after normalising whitespace and dropping data.table's column-type row, which only data.table ≥ 1.15 prints. Tables must match byte for byte.

## When a number changes

Don't overwrite it silently. In the folder README, add a dated note: what changed, why, and old → new. Then update the top-level README, and the post if it quotes the number. `imputation-hye/README.md` has an example.
