# FragPipe 24.0 × PXD028735 HYE (same 6 Orbitrap runs as ../maxquant-hye) — 2026-10-07

Tools (academic, non-commercial licences): FragPipe-24.0-linux, MSFragger-4.4.1, IonQuant-1.11.20, diaTracer-2.2.1
(diaTracer is unused for DDA but FragPipe 24 headless refuses to start without it). Java 21 Temurin in ~/.jdk.

**2026-11-23: re-run on verified raw files.** The first PRIDE downloads of A_Alpha_01 and A_Gamma_01 each contained a ~200-byte server error message (one damaged spectrum per file). FragPipe was re-run on verified copies with the same workflow, and its outputs came out identical (combined_protein/peptide/ion byte for byte), so every FragPipe number below is unchanged. Only `compare_mq_fp.out` changed, because it reads the re-run MaxQuant table: proteins in both 3,877 → 3,884, MaxQuant within ±0.5 98.9/92.2/80.9 → 98.9/91.5/82.6 %, FragPipe-only 2,905 → 2,898.

Pitfalls hit:
- FragPipe launches sub-tools with the `java` on PATH; system Java 8 → UnsupportedClassVersionError. Put Java 21 first on PATH.
- ThermoRawFileParser `-L` is --msLevel (not log level): `-L 1` wrote MS1-only mzML → MSFragger "Scans = 0".

Run: workflow = LFQ-MBR defaults (+db path, 16 threads); FASTA = same UniProt SP human/yeast/E. coli +
Philosopher contaminants + rev_ decoys (63,482). Workflow file: hye_lfq_mbr.workflow; manifest: hye.fp-manifest. The full run log is not included.
Score: Rscript score_fragpipe.R [workdir] [out dir] → fragpipe_scores.csv in the out dir, default the current folder (same rules as MaxQuant scoring; before 2026-11-02 it wrote to workdir/..).

RESULT (2026-10-07 15:49→15:59): ALL JOBS DONE IN 8.9 MINUTES (MSFragger 3.4, MSBooster 1.4, Percolator 1.2, IonQuant 1.6) · peak RSS 20 GB.
combined_protein.tsv 7,937 groups (human 5,069 · yeast 2,244 · E. coli 618). Scored MaxLFQ (≥2/3 per condition):
human/yeast/E. coli n 4,605/1,817/359, within ±0.5 94.7/88.3/73.5 %. Intensity: 84.1/72.5/44.7 %.
Same 3,884 proteins as MaxQuant LFQ: within ±0.5 MaxQuant 98.9/91.5/82.6 vs FragPipe 99.6/99.3/97.7 %;
centred on human median 99.4/96.6/82.6 vs 99.6/99.3/97.7. FragPipe-only scoreable 2,898 (86.7/78.1/59.5 %).
Note: Philosopher --contam entries carry no contam_ prefix here (0 flagged) — contaminants not removed for FragPipe.
Note on `compare_mq_fp.R` (the 3,884-protein comparison): species comes from the FragPipe `Protein` column alone, taking the first match in the order HUMAN, YEAST, ECOLI. `Indistinguishable Proteins` is not checked, and neither is the MaxQuant group's species. `score_fragpipe.R` instead requires one species across both columns. At most 6 FragPipe groups (7,937 − 7,931) lack a single species, so at most 6 of the 3,884 can differ.
