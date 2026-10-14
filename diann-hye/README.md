# DIA-NN 2.3.2 × PXD028735 HYE benchmark (2026-10-12/13)

Data: PRIDE PXD028735 (Van Puyvelde et al., Sci Data 2022, CC BY 4.0) — **SCIEX TripleTOF 6600+ SWATH** (99 variable windows). The dataset's Orbitrap DIA runs (AIF with staggered windows, GP-DIA) were not used.
Same design as the MaxQuant/FragPipe posts: Condition_{A,B}_Sample_{Alpha,Beta,Gamma}_01 → 3 vs 3. Truth log2(B/A): human 0, yeast −1, E. coli +2.
FASTA: same Swiss-Prot human+yeast+E. coli file as the FragPipe run (31,695 entries; fp-hye/db/hye_sp.fasta).
Raw files (30.3 GB, 6 .wiff + 6 .wiff.scan) are not in this repository; download them from PRIDE. report.pg_matrix.tsv and the derived matrices are included.

Download: pdl.py (10 ranged connections per file). PRIDE throttled to ~27 KB/s per connection and refused new connections for long stretches;
5 of 6 .wiff.scan files hit the 60-retry limit and were resumed with pdl_retry.py (2000 retries). ~8 h total.

Run (Windows DIA-NN 2.3.2 called from WSL, 24 threads, i9-14900K; DIA-NN ran natively on the Windows host with 128 GB RAM, WSL had 64 GB):
1. Library: --fasta-search --predictor --gen-spec-lib, trypsin K*/R*, 1 missed cleavage, 7–30 aa, charge 2–4, m/z 300–1800, fixed CAM, var Met-ox (max 1)
   → 6,370,129 precursors, 31,635 proteins, 11 min (libgen.log in C:\diann-hye).
2. Search: --lib predicted --reanalyse (MBR) --rt-profiling --qvalue 0.01 --matrices, mass accuracy automatic → 6 h 17 min.
   Each .wiff took ~17.5 min just to load (first pass ~45 min/file, second pass ~18 min/file).
   Per run, 2nd pass: ~100k precursors, ~8,200 protein groups at 1% FDR (report.stats.tsv).
3. No-MBR comparison: same command without --reanalyse plus --use-quant (reuses first-pass .quant; 171 s). DIA-NN warns that reusing .quant
   with automatic mass accuracy can give different results → treat as indicative only (run_nombr.log).

compare_three.R also needs the MaxQuant and FragPipe outputs (not included; regenerate with the other two folders).
Scoring: score_diann.R — same rules as MaxQuant/FragPipe (≥2 of 3 per condition, mean log2, within ±0.5 of truth). Output: score_diann.out, diann_scores.csv.
pg_from_parquet.py rebuilds a PG.MaxLFQ matrix from report.parquet with Q≤0.01 & PG.Q≤0.01 as a cross-check of report.pg_matrix.tsv.
Note: first-pass report has PG.MaxLFQ = 0 everywhere (no protein quant in the first pass) — hence the --use-quant rerun.
