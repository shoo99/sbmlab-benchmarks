# Biomarker data-leakage test (2026-10-13)
Data: Bader et al. 2020 Mol Syst Biol (PMC7266499, CC BY per Europe PMC) Dataset EV3 (log10, 1,484 PGs × 197 CSF samples) + EV2 (labels, cohorts),
downloaded via Europe PMC supplementaryFiles (PMC's own links return a bot page). Raw xlsx not committed.
prep.py → X_log10.csv / labels.csv (not committed; regenerate). leak.py K → leak_k{5,20,100}.csv (50 repeats; shuffled = 50 permutations).
external.py → external.out / external_shuffled.csv (train Berlin+Sweden, test Magdeburg+Kiel). chart.py → sbmlab figure.
post_code.py = the code block in the post, run as written: (197, 1484) leaky 0.72 (max 0.83) | proper 0.50 (max 0.67).
sklearn 1.9.1, pandas 3.0.6, numpy 2.5.3.
