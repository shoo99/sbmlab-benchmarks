# Bioinformatics database API check (2026-10-13, Seoul, ~23:00 KST)
One protein (human SOD1, UniProt P00441, NCBI Gene 6647) through 13 public APIs.
probe.py N: each endpoint N times, 1 s apart. Round 1 -> probe_run1.csv/.out; round 2 -> probe.csv/.out (probe.json with the raw response bodies is not included: some sources, e.g. KEGG, do not allow redistribution; rerun probe.py to regenerate it).
extract.py: what each DB returned (from probe.json). versions.json: Ensembl 116, STRING 12.5, Reactome 97, KEGG pathway last update (fetched 23:20).
one_protein.py = the code block in the sbmlab post, run as written (output in the post; 400/404 raise immediately, verified).
GEO 526 = 70 Series + 456 Samples (gse/gsm[ETYP]).
