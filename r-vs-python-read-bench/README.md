# R vs Python read benchmark (2026-10-10)
bench.R / bench.py on MaxQuant 2.8.1 evidence.txt (HYE, 146 MB, 271,777 rows). bench.out = results.
post_snippets_check.py: the sbmlab post's nested-CV and polars snippets run on synthetic data (CV accuracy 0.783).

## Note (2026-11-02): base R did more work than the others
- The base R timing in `bench.out` (5.80 s) is `read.delim` as written in the post. It parses and type-converts every column of evidence.txt, while data.table (`select=`), pandas (`usecols`) and polars (`.select`) read only the 4 columns the summary needs. So the ~100× gap to data.table overstates the difference in reading speed. Base R's `aggregate()` also drops rows with NA Intensity.
- `bench.R` now also times base R reading the same 4 columns (`colClasses = "NULL"` for the rest, NA rows kept, mean and n per protein). It also prints that variant's group count, which should equal data.table's. `bench.out` predates this and was not re-run, because evidence.txt is not in the repository. As a rough check only: on a synthetic 60,000-row, 40-column file, the 4-column read took about a quarter of the all-column time.
- Threads: fread used 16 threads and polars 32 (both printed in `bench.out`); pandas `read_csv` is single-threaded.
