# STRING enrichment background test (2026-10-13)
**2026-11-23: re-run on the re-searched DIA-NN matrix.** `../diann-hye/report.pg_matrix.tsv` was replaced when the SWATH runs were searched again on verified raw files (the first downloads had server error messages written into 5 of 6 .wiff.scan files; see diann-hye). The detected set changed slightly (5,661 → 5,640 human groups), so this folder was re-run. With the old matrix every script reproduced its bundled output byte for byte (counts.csv, planted_terms.out, post_code.R 115 / 0). Median terms for random lists, genome background: 103 → 98.5 (range 63–168 → 51 to 186); detected background: 0 in every list, as before.

Detected set: human protein groups in ../diann-hye/report.pg_matrix.tsv (5,640 groups → 5,660 STRING v12.0 ids).
run_string_bg.R: 20 random draws of 300 detected proteins + 5 planted lists (15 Reactome cholesterol-biosynthesis + 285 random),
each enriched via STRINGdb 2.14.3 get_enrichment() with default (genome) vs set_background(detected). Outputs counts.csv,
null_terms_genome_background.csv, planted.csv, run.log. planted_terms.R → planted_terms.out (terms + rank 15 of 111 for the first sterol term; 21 of 110 before the re-run).
post_code.R = the code block in the sbmlab post, run as written (seed 1: genome 118, detected 0, matching counts.csv; 115 before the re-run).
chart.py → sbmlab public/images/blog/string-background-null.png.
