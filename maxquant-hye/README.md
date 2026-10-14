# MaxQuant 2.8.1 × PXD028735 HYE benchmark (2026-10-07)

Data: PRIDE PXD028735 (Van Puyvelde et al., Sci Data 2022, CC BY 4.0) — Orbitrap QE HF-X DDA, 120 min.
Condition A = Human 65 / Yeast 30 / E. coli 5 ; B = 65 / 15 / 20 → log2(B/A): 0 / −1 / +2.
Runs used: Condition_{A,B}_Sample_{Alpha,Beta,Gamma}_01 (one injection per preparation → 3 vs 3).
FASTA: UniProt Swiss-Prot human 9606 (20,431) · yeast 559292 (6,733) · E. coli 83333 (4,531), fetched 2026-10-06.

Run: Linux (WSL Ubuntu 24.04), .NET 8.0.31, MaxQuantCmd 2.8.1.0, 16 threads, LFQ (min ratio count 2) + MBR + iBAQ.
Wall time 1:01:34, peak RSS ~11 GB (main process). First attempt (28 threads) died with WSL, most likely because the Windows C: drive filled.
Raw files (21.7 GB) and the full MaxQuant output are not in this repository; download the raw files from PRIDE.

Files: mqpar.xml (as run) · score_quant.R → quant_scores.csv · score_mbr.R → mbr_scores.csv ·
summary.txt · parameters.txt · #runningTimes.txt · pdl.py (PRIDE ranged downloader; PRIDE caps ~0.35 MB/s per
connection) · watchdog.sh (stops MaxQuant before WSL disk / C: / RAM run out).

Findings (written up in the sbmlab.com post linked from the top-level README):
- MaxQuant 2.8 names the decoy column `Decoy` (was `Reverse`): base R `pg[pg$Reverse != "+", ]` returns 0 rows silently.
- `MaxQuantCmd --create` truncated experiment names (Condition_A_Sample_Alpha_01 → A_Sample_Alph).
- Filters: Decoy 68 · contaminant 59 · only-by-site 44 → 5,628 groups (human 3,780 · yeast 1,457 · E. coli 389 · mixed 2).
- iBAQ between-sample ratios are identical to Intensity ratios.
- Within ±0.5 of truth (human/yeast/E. coli): Intensity 86/74/43 % vs LFQ 99/92/81 %; LFQ quantifies fewer (2,896/907/132 vs 3,603/1,250/243). LFQ human median +0.14 (normalisation assumption violated by design).
- MBR filled 20,830 of 210,449 peptide cells (9.9 %); accuracy drops with matching, most for low-abundance E. coli (88 → 77 → 62 %).

MBR off rerun (2026-10-07 12:38→13:30, 52.0 min): LFQ scoreable human/yeast/E. coli 2,346/656/77 (on 2,896/907/132);
same-protein within ±0.5 on vs off: 99.5/99.1 · 97.7/95.0 · 94.8/62.3 %; MBR-only proteins 860 (96/78/62 %).
Intensity within ±0.5 on/off: 86/77 · 74/59 · 43/20 %. Script compare_mbr_onoff.R → mbr_onoff_*.csv
