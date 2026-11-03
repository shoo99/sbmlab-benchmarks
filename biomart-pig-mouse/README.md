# BioMart pig->mouse 1:1 orthologs (2026-10-10)
biomart_query.py → biomart_query.out: release 111 (jan2024 archive) 16,510 one2one; 116 (jun2026 archive) 16,529; shared 16,004.
r_fallback.R: plain-R martservice query (the post's R code, plus a completeness check added 2026-11-02: `completionStamp="1"` and a test for the closing `[success]` line, which biomart_query.py now also requires; without it a truncated or failed BioMart answer is parsed silently); it returned 26,676 rows / 16,510 one2one when run (output not saved as a file; same release-111 count as biomart_query.out).
biomart_counts.R: biomaRt 2.58.2 attempt — fails (invalid 'file' argument listing archives; www.ensembl.org/biomart 308 → jun2026 archive).
pybiomart: host must be http:// (https:// → name resolution error); use_attr_names=True works.
Mixing uniprotswissprot with homolog attributes → 'Attributes from multiple attribute pages are not allowed'.
