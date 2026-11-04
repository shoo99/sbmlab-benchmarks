# Metagenomic taxonomic profilers on ZymoBIOMICS mock communities (in progress)

Tracking issue: [#15](https://github.com/shoo99/sbmlab-benchmarks/issues/15). No post yet.

**Status: running.** The scripts and the truth table are in place. Kraken2 + Bracken with two capped databases, standard-8 (no fungi) and PlusPF-8 (with fungi), are being run on two mock communities:

- ZymoBIOMICS D6300, even: ERR2984773
- ZymoBIOMICS D6310, log distribution: ERR2935805

MetaPhlAn 4, mOTUs and sylph are workstation steps, because their databases are not reachable from the cloud environment.

**Truth** (`truth.tsv`): the theoretical genomic-DNA composition from the Zymo product sheets.

- D6300: 8 bacteria at 12% each and 2 yeasts at 2% each.
- D6310: 89.1% down to 0.000089%.

**Scripts**:

- `download.sh` fetches the inputs.
- `run.sh kraken` runs Kraken2 and Bracken, followed by the workstation steps `metaphlan`, `motus` and `sylph`.
- `score.py` writes the scores, and `figure.py` draws the figure.
