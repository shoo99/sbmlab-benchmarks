# DE pipeline post check (2026-10-10)
post_code.R: every R block of sbmlab differential-expression-proteomics-complete-r-pipeline, run on the HYE MaxQuant 2.8.1 LFQ table (test-harness lines marked). Symbol-based enrichGO block skipped (no symbols); UniProt block tested with human IDs.

**2026-11-23: re-run on verified raw files.** The first PRIDE downloads of some raw files contained short server error messages, so MaxQuant and DIA-NN were re-run on verified copies with the same settings. `post_code.R` keeps no saved output here and still reads the first MaxQuant run (~/mq-hye/results_mbr_on/txt); the numbers in `diann270/` are from the re-run (see its README).
