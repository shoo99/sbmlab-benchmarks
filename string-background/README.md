# STRING enrichment background test (2026-10-13)
Detected set: human protein groups in ../diann-hye/report.pg_matrix.tsv (5,661 groups → 5,686 STRING v12.0 ids).
run_string_bg.R: 20 random draws of 300 detected proteins + 5 planted lists (15 Reactome cholesterol-biosynthesis + 285 random),
each enriched via STRINGdb 2.14.3 get_enrichment() with default (genome) vs set_background(detected). Outputs counts.csv,
null_terms_genome_background.csv, planted.csv, run.log. planted_terms.R → planted_terms.out (terms + rank 21 of first sterol term).
post_code.R = the code block in the sbmlab post, run as written (seed 1: genome 115, detected 0 — matches counts.csv).
chart.py → sbmlab public/images/blog/string-background-null.png.
