# R vs Python read benchmark (2026-10-10)
bench.R / bench.py on MaxQuant 2.8.1 evidence.txt (HYE, 146 MB, 271,777 rows). bench.out = results.
post_snippets_check.py: the sbmlab post's nested-CV and polars snippets run on synthetic data (CV accuracy 0.783).

## Note (2026-11-02): base R did more work than the others
- The base R timing in `bench.out` (5.80 s) is `read.delim` as written in the post. It parses and type-converts every column of evidence.txt, while data.table (`select=`), pandas (`usecols`) and polars (`.select`) read only the 4 columns the summary needs. So the ~100× gap to data.table overstates the difference in reading speed. Base R's `aggregate()` also drops rows with NA Intensity.
- `bench.R` now also times base R reading the same 4 columns (`colClasses = "NULL"` for the rest, NA rows kept, mean and n per protein). It also prints that variant's group count, which should equal data.table's. `bench.out` was re-run on 2026-11-03 (same evidence.txt, same machine, pandas 3.1.0 and polars 2.1.0): data.table 0.07 s, base R all columns 8.22 s, base R 4 columns 2.25 s, pandas 0.37 s, polars 0.15 s. Reading the same 4 columns, base R is about 30 times slower than data.table, not about 100 times. The 2026-10-10 run (data.table 0.04 s, base R 5.80 s, pandas 3.0.6 0.51 s, polars 2.0.0 0.16 s) is in the git history; compare timings only within one run.
- Threads: fread used 16 threads and polars 32 (both printed in `bench.out`); pandas `read_csv` is single-threaded.
