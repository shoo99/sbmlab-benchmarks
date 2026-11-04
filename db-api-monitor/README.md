# Long-term availability of bioinformatics database APIs (weekly check, since 2026-11)

Tracking issue: [#9](https://github.com/shoo99/sbmlab-benchmarks/issues/9). No post yet.

How reliable are the public APIs that people script against, measured over months rather than on one evening? [`bioinfo-db-check`](../bioinfo-db-check) followed one protein (human SOD1) through 13 APIs once, on 2026-10-13. This folder repeats a light version of that check every week from GitHub Actions and keeps every result.

**Status: the weekly check starts after this folder is merged.** `checks.csv` holds only its header until the first run. Until then, `scores.csv` holds only the 2026-10-13 evening measurement, which stays as the local reference.

## What is measured

The 13 endpoints of `bioinfo-db-check/probe.py`:

- UniProt, Ensembl, NCBI Gene, ClinVar
- RCSB PDB, AlphaFold DB, InterPro, QuickGO
- STRING, Reactome, KEGG, GEO, PRIDE

Each is called **once per run**, 1 s apart. `monitor.py` appends one row per API to `checks.csv`:

| Column | Content |
|---|---|
| `date_utc` | start of the run (UTC) |
| `runner` | where it ran (`github-actions ubuntu-24.04`) |
| `db` | database |
| `status` | HTTP status, or `error:timeout` / `error:dns` / `error:connection` / ... |
| `seconds` | wall-clock time of the main request |
| `release` | release or version the API reports |
| `example_ok` | 1 if the response still holds the documented example, else 0 |

Where the release comes from:

| Database | Release source |
|---|---|
| UniProt | `X-UniProt-Release` header of the main response |
| Ensembl | `/info/data` |
| AlphaFold DB | `latestVersion` in the main response |
| InterPro | the API root |
| STRING | `/api/json/version` |
| Reactome | `/data/database/version` |
| KEGG | `/info/pathway` |

Each of the five version requests is one small extra call. The other six APIs report no release this cheaply.

**The documented example** is what `bioinfo-db-check` saw. UniProt must return P00441, Ensembl ENSG00000142168, NCBI Gene 6647, AlphaFold `AF-P00441-…`, and KEGG at least one `path:hsa` link. ClinVar, PDB, InterPro, QuickGO, STRING, Reactome, GEO and PRIDE must each return at least one hit. A response that is HTTP 200 but no longer holds the example counts as a broken example.

## Rules

- **No response bodies are stored**, only the columns above. Some sources, KEGG for example, do not allow redistribution.
- **Light on the APIs**: one call per API per week, plus the five version requests, 1 s apart.
- **User-Agent** names the project and links the repository, without a personal contact: `sbmlab-benchmarks-db-api-monitor/1.0 (+https://github.com/shoo99/sbmlab-benchmarks)`.
- **Latency is not what a user in Korea sees.** GitHub's runners sit in a few cloud regions, mostly in the US. The 2026-10-13 measurement from a home network in Seoul stays as the local reference, in its own rows of `scores.csv`.

## Scores

`summarize.py` uses the standard library only. It writes these metrics per database to `scores.csv`:

| Metric | Meaning |
|---|---|
| `uptime_share` | runs with HTTP 200 |
| `example_ok_share` | runs whose response still held the example |
| `median_seconds`, `p95_seconds` | latency of the successful runs (p95 by nearest rank) |
| `release_changes` | how often the reported release changed |

Its printout, `summarize.out`, also lists the dates of release changes, of broken examples and of failed calls.

## How it runs

- **Workflow**: [`db-api-monitor.yml`](../.github/workflows/db-api-monitor.yml) runs on Wednesdays at 02:29 UTC (11:29 KST), and on demand (Actions → db-api-monitor → Run workflow).
- **What it does**:
  1. `monitor.py` checks the APIs.
  2. The new rows are appended to `checks.csv`.
  3. `summarize.py` and `docs/make_summaries.py` refresh `scores.csv`, `summarize.out` and the generated index pages.
  4. `github-actions[bot]` commits the result to `main` (`contents: write`). If `main` moved in the meantime, it starts again from the new `main`, so generated files are never merged by hand.
- **If `main` is protected** against direct pushes, this step fails. Allow `github-actions[bot]` to push, or change the step to open a pull request.
- **Run it by hand**: `python3 monitor.py` appends to `checks.csv` (set `MONITOR_RUNNER` to say where it ran), then run `python3 summarize.py > summarize.out`.
- **CI**:
  - re-runs `summarize.py` and compares `scores.csv` and `summarize.out` with the bundled files
  - runs `test/test_monitor.py`. It points the 13 endpoints at a local HTTP server that answers like the real APIs, then checks 20 known answers: statuses, every release parser, example checks, an unreachable host, that no body is stored, and the summary's uptime and release changes.

## Results

Pending the first weekly runs. The local reference from `bioinfo-db-check` (2026-10-13, about 23:00 KST, 2 rounds × 5 calls) is in `summarize.out`. For example:

| Database | Uptime | Median latency |
|---|---|---|
| Ensembl | 6 of 10 calls | 8.1 s |
| KEGG | 9 of 10 | 0.28 s |
| PRIDE | 10 of 10 | 7.3 s |

## Files

| File | Content |
|---|---|
| `monitor.py` | one check of the 13 APIs; appends to `checks.csv` |
| `checks.csv` | every weekly result (status, seconds, release, example check) |
| `summarize.py`, `summarize.out`, `scores.csv` | per-database uptime, latency and release changes, plus the local reference |
| `test/test_monitor.py` | offline known-answer test (CI) |
