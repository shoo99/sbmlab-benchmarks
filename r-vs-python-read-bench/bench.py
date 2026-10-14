# Same task in pandas and polars; median of 5 runs.
import sys, time, statistics, pandas as pd, polars as pl
f = sys.argv[1]; cols = ["Leading razor proteins", "PEP", "Intensity", "Decoy"]
def t(fn): r = []; [r.append((lambda s: (fn(), time.perf_counter() - s)[1])(time.perf_counter())) for _ in range(5)]; return statistics.median(r)
def pd_run():
    d = pd.read_csv(f, sep="\t", usecols=cols, quoting=3, low_memory=False)
    d = d[(d["Decoy"] != "+") & (d["PEP"] < 0.01)]
    return d.groupby("Leading razor proteins")["Intensity"].agg(["mean", "size"])
def pl_run():
    return (pl.scan_csv(f, separator="\t", quote_char=None, infer_schema_length=10000)
              .select(cols).filter((pl.col("Decoy").fill_null("") != "+") & (pl.col("PEP") < 0.01))
              .group_by("Leading razor proteins").agg(pl.col("Intensity").mean().alias("mean_int"), pl.len().alias("n")).collect())
print("groups pandas", len(pd_run()), "polars", pl_run().height)
print(f"pandas {pd.__version__} read_csv: median {t(pd_run):.2f} s")
print(f"polars {pl.__version__} scan_csv (threads {pl.thread_pool_size()}): median {t(pl_run):.2f} s")
