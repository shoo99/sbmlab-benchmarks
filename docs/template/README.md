# <Title: the question this folder answers> (<YYYY-MM-DD>)

Post: [<post title>](https://sbmlab.com/blog/<slug>)

<One or two sentences: what was measured and the headline result, with numbers.>

## Data

- Source: <repository and accession, or paper and supplementary table>, downloaded <date>.
- Licence: <licence of the data>. <Whether the raw files are in this folder; usually not.>
- Download: <command, or the steps to get the files>.

## Truth

<What the right answer is and how it is known: known mixing ratios, spike-ins, a curated truth set, a simulation with planted effects, shuffled labels whose expected result is chance, or "none" when the folder measures agreement or speed instead. State any centring or normalisation applied before scoring.>

## Tools

| Tool | Version | Settings that matter |
|---|---|---|
| <tool> | <x.y.z> | <non-default parameters> |

Machine: <CPU, RAM, OS / WSL>. R <version>, Python <version>.

## How to run

```bash
./run.sh <input dir>
```

| Script | What it does | Output |
|---|---|---|
| `<script>` | <one line> | `<script>.out`, `<table>.csv` |

Outputs are written next to the scripts. Inputs are passed as arguments, so nothing in the scripts points to one machine.

## Results

| <Metric> | <Tool A> | <Tool B> |
|---|---|---|
| <...> | <...> | <...> |

Every number above is in `<file>` (line or column <...>). New folders also write `scores.csv` in the standard long format (see CONTRIBUTING.md).

## Traps

- <What went wrong or surprised you, and how you found out.>

## Not included

- <Raw data, files that may not be redistributed, intermediates over ~10 MB (link to the Zenodo record instead).>

## Corrections

<Dated notes if a number changes after the post: what changed, why, old → new. Leave out until needed.>
