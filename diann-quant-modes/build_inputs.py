"""Build protein-level matrices for DIA-NN quantification modes from the same .quant files (PXD028735 HYE, SWATH).
Usage: python build_inputs.py [first|nombr|mbr|nombr_runq]
  first (default): first exploratory reruns (--use-quant, no MBR, mass accuracy auto):
    out_frtest = QuantUMS high-precision (default), out_highacc = --high-acc, out_direct = --direct-quant (legacy); writes here.
  nombr / mbr: clean rerun (clean_rerun.sh: mass accuracy fixed 16/13 ppm, window 9) with MBR off (--use-quant)
    or on (--reanalyse): c_<set>_quantums / c_<set>_highacc / c_<set>_direct; writes to clean/<set>/.
  Modes whose run has not finished (no report.pg_matrix.tsv) are skipped.
Writes:
  pg_<mode>.tsv              DIA-NN's own report.pg_matrix.tsv, renamed columns A_Alpha..B_Gamma
  pg_precsum_<mode>.tsv      sum of Precursor.Normalised per protein group
  pg_precsum6_<mode>.tsv     the same sum restricted to precursors quantified in all six runs
  Precursor filter: first = Q.Value & PG.Q.Value <= 0.01; clean sets = the filter that reproduces DIA-NN's own
  report.pr_matrix.tsv exactly (checked by cell count): nombr = Global.Q.Value & Global.PG.Q.Value <= 0.01,
  mbr = Lib.Q.Value & Lib.PG.Q.Value <= 0.01. nombr_runq = the MBR-off runs with the run-specific filter
  (Q.Value & PG.Q.Value <= 0.01) instead, kept only to show how much the filter alone moves the fragment routes.
  fragments_long.tsv         fragment-level apex quantities (Fr.N.Quantity, x normalisation factor) for iq / directLFQ
"""
import sys, os
import numpy as np, pandas as pd, pyarrow.parquet as pq
BASE = "/mnt/c/diann-hye"
SET = sys.argv[1] if len(sys.argv) > 1 else "first"
RUNS = {"A_Sample_Alpha": "A_Alpha", "A_Sample_Beta": "A_Beta", "A_Sample_Gamma": "A_Gamma",
        "B_Sample_Alpha": "B_Alpha", "B_Sample_Beta": "B_Beta", "B_Sample_Gamma": "B_Gamma"}
def short(run):
    for k, v in RUNS.items():
        if k in run: return v
    raise ValueError(run)
if SET == "first":
    MODES = {"quantums_precision": "out_frtest", "quantums_accuracy": "out_highacc", "legacy_direct": "out_direct"}
else:
    run = SET.split("_")[0]   # nombr_runq = the MBR-off runs with the run-specific filter (to show what the filter does)
    MODES = {"quantums_precision": f"c_{run}_quantums", "quantums_accuracy": f"c_{run}_highacc", "legacy_direct": f"c_{run}_direct"}
    os.makedirs(f"clean/{SET}", exist_ok=True); os.chdir(f"clean/{SET}")
QC = {"first": ("Q.Value", "PG.Q.Value"), "nombr": ("Global.Q.Value", "Global.PG.Q.Value"), "mbr": ("Lib.Q.Value", "Lib.PG.Q.Value"),
      "nombr_runq": ("Q.Value", "PG.Q.Value")}[SET]
def keep(r): return r[(r[QC[0]] <= 0.01) & (r[QC[1]] <= 0.01) & (r["Decoy"] == 0)]
for mode, d in MODES.items():
    if not os.path.exists(f"{BASE}/{d}/report.pg_matrix.tsv"):
        print(mode, "skipped (not finished)"); continue
    pg = pd.read_csv(f"{BASE}/{d}/report.pg_matrix.tsv", sep="\t")
    cols = [c for c in pg.columns if c.endswith(".wiff")]
    pg = pg.rename(columns={c: short(c) for c in cols})
    pg[["Protein.Group", "Protein.Names"] + list(RUNS.values())].to_csv(f"pg_{mode}.tsv", sep="\t", index=False)
    r = pq.read_table(f"{BASE}/{d}/report.parquet", columns=["Run", "Protein.Group", "Protein.Names", "Precursor.Id", "Q.Value", "PG.Q.Value",
                                                             "Precursor.Normalised", "Precursor.Quantity", "Decoy"] + [c for c in QC if c not in ("Q.Value", "PG.Q.Value")]).to_pandas()
    r = keep(r)
    pm = pd.read_csv(f"{BASE}/{d}/report.pr_matrix.tsv", sep="\t")
    print(mode, "kept precursor rows", len(r), "| non-missing cells in DIA-NN pr_matrix", int(pm[[c for c in pm.columns if c.endswith(".wiff")]].notna().sum().sum()))
    r["sample"] = r["Run"].map(short)
    s = r.pivot_table(index=["Protein.Group", "Protein.Names"], columns="sample", values="Precursor.Normalised", aggfunc="sum").reset_index()
    s.to_csv(f"pg_precsum_{mode}.tsv", sep="\t", index=False)
    # same sum, but only precursors quantified in all six runs (so a protein's sum uses the same precursors in A and B)
    w = r.pivot_table(index=["Protein.Group", "Protein.Names", "Precursor.Id"], columns="sample", values="Precursor.Normalised", aggfunc="sum")
    w = w[w.notna().all(axis=1)]
    w.groupby(level=[0, 1]).sum().reset_index().to_csv(f"pg_precsum6_{mode}.tsv", sep="\t", index=False)
    print(mode, "pg_matrix rows", len(pg), "| precursor rows", len(r), "| precursor-sum proteins", len(s))
# fragment-level from the default run (fragment apex quantities do not depend on the protein-level mode)
cols = ["Run", "Protein.Group", "Protein.Names", "Precursor.Id", "Q.Value", "PG.Q.Value", "Decoy", "Precursor.Normalised", "Precursor.Quantity"]
cols += [c for c in QC if c not in cols]
cols += [f"Fr.{i}.Quantity" for i in range(12)] + [f"Fr.{i}.Id" for i in range(12)]
r = pq.read_table(f"{BASE}/{MODES['quantums_precision']}/report.parquet", columns=cols).to_pandas()
r = keep(r).copy()
r["sample"] = r["Run"].map(short)
nf = np.where(r["Precursor.Quantity"] > 0, r["Precursor.Normalised"] / r["Precursor.Quantity"], np.nan)  # run normalisation factor
parts = []
for i in range(12):
    p = r[["sample", "Protein.Group", "Protein.Names", "Precursor.Id"]].copy()
    p["q"] = r[f"Fr.{i}.Quantity"].to_numpy() * nf
    p["ion"] = r["Precursor.Id"] + "|" + r[f"Fr.{i}.Id"].astype(str)
    parts.append(p)
long = pd.concat(parts)
long = long[long["q"] > 0]
long[["sample", "Protein.Group", "Protein.Names", "ion", "q"]].to_csv("fragments_long.tsv", sep="\t", index=False)
print("fragment rows", len(long), "| proteins", long["Protein.Group"].nunique(), "| ions", long["ion"].nunique())
