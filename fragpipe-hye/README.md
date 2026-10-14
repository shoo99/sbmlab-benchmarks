# FragPipe 24.0 × PXD028735 HYE (same 6 Orbitrap runs as ../maxquant-hye) — 2026-10-07

Tools (academic, non-commercial licences): FragPipe-24.0-linux, MSFragger-4.4.1, IonQuant-1.11.20, diaTracer-2.2.1
(diaTracer is unused for DDA but FragPipe 24 headless refuses to start without it). Java 21 Temurin in ~/.jdk.

Pitfalls hit:
- FragPipe launches sub-tools with the `java` on PATH; system Java 8 → UnsupportedClassVersionError. Put Java 21 first on PATH.
- ThermoRawFileParser `-L` is --msLevel (not log level): `-L 1` wrote MS1-only mzML → MSFragger "Scans = 0".

Run: workflow = LFQ-MBR defaults (+db path, 16 threads); FASTA = same UniProt SP human/yeast/E. coli +
Philosopher contaminants + rev_ decoys (63,482). Workflow file: hye_lfq_mbr.workflow; manifest: hye.fp-manifest. The full run log is not included.
Score: Rscript score_fragpipe.R [workdir] → fragpipe_scores.csv (same rules as MaxQuant scoring).

RESULT (2026-10-07 15:49→15:59): ALL JOBS DONE IN 8.9 MINUTES (MSFragger 3.4, MSBooster 1.4, Percolator 1.2, IonQuant 1.6) · peak RSS 20 GB.
combined_protein.tsv 7,937 groups (human 5,069 · yeast 2,244 · E. coli 618). Scored MaxLFQ (≥2/3 per condition):
human/yeast/E. coli n 4,605/1,817/359, within ±0.5 94.7/88.3/73.5 %. Intensity: 84.1/72.5/44.7 %.
Same 3,877 proteins as MaxQuant LFQ: within ±0.5 MaxQuant 98.9/92.2/80.9 vs FragPipe 99.7/99.2/97.7 %;
centred on human median 99.4/96.9/82.4 vs 99.7/99.2/97.7. FragPipe-only scoreable 2,905 (86.6/78.3/59.6 %).
Note: Philosopher --contam entries carry no contam_ prefix here (0 flagged) — contaminants not removed for FragPipe.
