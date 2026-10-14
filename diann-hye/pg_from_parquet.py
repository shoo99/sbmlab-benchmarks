"""Build a protein-group matrix (PG.MaxLFQ) from a DIA-NN report parquet, one filter for both passes:
precursor Q.Value <= 0.01 and run-level PG.Q.Value <= 0.01. Output columns match report.pg_matrix.tsv."""
import sys, polars as pl
src, out = sys.argv[1], sys.argv[2]
d = pl.read_parquet(src, columns=["Run", "Protein.Group", "Protein.Names", "PG.MaxLFQ", "Q.Value", "PG.Q.Value"])
d = d.filter((pl.col("Q.Value") <= 0.01) & (pl.col("PG.Q.Value") <= 0.01))
d = d.group_by(["Protein.Group", "Protein.Names", "Run"]).agg(pl.col("PG.MaxLFQ").first())
m = d.pivot(on="Run", index=["Protein.Group", "Protein.Names"], values="PG.MaxLFQ")
m.write_csv(out, separator="\t")
print(src, "->", out, m.shape)
