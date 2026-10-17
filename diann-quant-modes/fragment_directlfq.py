"""Fragment-level directLFQ (Ammar et al. 2023) on the same fragment table."""
import sys, os, pandas as pd, directlfq.lfq_manager as lm
if len(sys.argv) > 1: os.chdir(sys.argv[1])   # e.g. clean/nombr (default: this folder = first reruns)
d = pd.read_csv("fragments_long.tsv", sep="\t")
w = d.pivot_table(index=["Protein.Group", "ion"], columns="sample", values="q", aggfunc="max").reset_index()
w = w.rename(columns={"Protein.Group": "protein"})
w.to_csv("fragments_wide.aq_reformat.tsv", sep="\t", index=False)
lm.run_lfq("fragments_wide.aq_reformat.tsv", num_cores=8)
