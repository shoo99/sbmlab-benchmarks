# MaxQuant 2.8.1 × PXD028735 HYE benchmark (2026-10-07)

Data: PRIDE PXD028735 (Van Puyvelde et al., Sci Data 2022, CC BY 4.0) — Orbitrap QE HF-X DDA, 120 min.
Condition A = Human 65 / Yeast 30 / E. coli 5 ; B = 65 / 15 / 20 → log2(B/A): 0 / −1 / +2.
Runs used: Condition_{A,B}_Sample_{Alpha,Beta,Gamma}_01 (one injection per preparation → 3 vs 3).
FASTA: UniProt Swiss-Prot human 9606 (20,431) · yeast 559292 (6,733) · E. coli 83333 (4,531), fetched 2026-10-06.

**2026-11-23: re-run on verified raw files.** The first PRIDE downloads of A_Alpha_01 and A_Gamma_01 each contained a ~200-byte server error message (one damaged spectrum per file). Both searches (MBR on and off) were repeated on verified copies with the same mqpar settings, and every number below is from the re-run. Old → new: LFQ within ±0.5 human/yeast/E. coli 98.8/92.3/81.1 → 98.8/91.6/82.7 %; LFQ scoreable 2,896/907/132 → 2,896/912/133 (MBR on) and 2,346/656/77 → 2,365/664/75 (MBR off). The old MBR-off run had different per-run peptide counts from the MBR-on run; in the re-run they match exactly (`summary_mbr_off.txt` = `summary.txt`).

Run: Linux (WSL Ubuntu 24.04), .NET 8.0.31, MaxQuantCmd 2.8.1.0, 16 threads, LFQ (min ratio count 2) + MBR + iBAQ.
Wall time 1:01:34, peak RSS ~11 GB (main process). The re-run (2026-11-19) took 8:49:11 (`#runningTimes.txt`), with the MBR-off run going at the same time and the raw files on a Windows drive mounted in WSL (/mnt/g), so its time is not comparable. First attempt (28 threads) died with WSL, most likely because the Windows C: drive filled.
Raw files (21.7 GB) and the full MaxQuant output are not in this repository; download the raw files from PRIDE.

Files: mqpar.xml (as run) · score_quant.R → quant_scores.csv · score_mbr.R → mbr_scores.csv ·
summary.txt · parameters.txt · #runningTimes.txt · pdl.py (PRIDE ranged downloader; PRIDE caps ~0.35 MB/s per
connection) · watchdog.sh (stops MaxQuant before WSL disk / C: / RAM run out).

Findings (written up in the sbmlab.com post linked from the top-level README):
- MaxQuant 2.8 names the decoy column `Decoy` (was `Reverse`): base R `pg[pg$Reverse != "+", ]` returns 0 rows silently.
- `MaxQuantCmd --create` truncated experiment names (Condition_A_Sample_Alpha_01 → A_Sample_Alph).
- Filters: Decoy 68 · contaminant 58 · only-by-site 44 → 5,628 groups (human 3,783 · yeast 1,456 · E. coli 387 · mixed 2).
- iBAQ between-sample ratios are identical to Intensity ratios.
- Within ±0.5 of truth (human/yeast/E. coli): Intensity 86/74/44 % vs LFQ 99/92/83 %; LFQ quantifies fewer (2,896/912/133 vs 3,601/1,253/243). LFQ human median +0.14 (normalisation assumption violated by design).
- MBR filled 20,814 of 210,627 peptide cells (9.9 %); accuracy drops with matching, most for low-abundance E. coli (88 → 76 → 62 %).

MBR off rerun (first run 2026-10-07 12:38→13:30, 52.0 min; re-run 2026-11-19, 467.2 min alongside MBR on): LFQ scoreable human/yeast/E. coli 2,365/664/75 (on 2,896/912/133);
same-protein within ±0.5 on vs off: 99.5/99.2 · 97.7/94.3 · 94.7/64.0 %; MBR-only proteins 837 (96/75/67 %).
Intensity within ±0.5 on/off: 86/77 · 74/59 · 44/22 %. Script compare_mbr_onoff.R → mbr_onoff_*.csv

## Running the scripts

- `Rscript score_quant.R [txt dir] [out dir]`, `Rscript score_mbr.R [txt dir] [out dir] [fasta dir]` and `Rscript compare_mbr_onoff.R [MBR-on proteinGroups.txt] [MBR-off proteinGroups.txt] [out dir]`. Outputs go to the out dir, by default the current folder. Before 2026-11-02 the scripts always wrote to `/home/sysoft/mq-hye/analysis/`, even when given another input folder.
- `compare_mbr_onoff.R` now also saves the two tables it used to only print: `mbr_onoff_same_proteins.csv` (the "same-protein" line above) and `mbr_onoff_only_on_<quant>.csv` (MBR-only proteins). Both were written on 2026-11-03 by re-running the script on the raw MaxQuant output; the two existing `mbr_onoff_<quant>.csv` files came out byte for byte identical, and the same-protein values match the printout above. `preprint-hye/figures.py` uses them for Fig 2b.
- Peptide species in `score_mbr.R` comes from the `Leading razor protein` only. A peptide shared between species is scored against its razor protein's species, not dropped. The protein-level scripts instead require all proteins in a group to come from one species. How many peptides map to more than one species was not counted.
