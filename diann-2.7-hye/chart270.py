"""Figure for the DIA-NN 2.7.0 vs 2.3.2 post: (left) share within +/-0.5 of the true log2 ratio by protein set,
(right) E. coli log2(B/A) against intensity in condition A, shared vs 2.7.0-only proteins."""
import sys, csv, statistics as st, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
rows = list(csv.DictReader(open("fc_both_versions.csv")))
for r in rows: r["fc"] = float(r["fc"]); r["intA"] = float(r["intA"])
T = {"Human": 0, "Yeast": -1, "Ecoli": 2}
def share(sel): return 100 * sum(abs(r["fc"] - T[r["species"]]) <= 0.5 for r in sel) / len(sel), len(sel)
groups = [("2.3.2, all", lambda r: r["version"] == "2.3.2"), ("2.7.0, all", lambda r: r["version"] == "2.7.0"),
          ("2.3.2, shared", lambda r: r["version"] == "2.3.2" and r["set"] == "both"),
          ("2.7.0, shared", lambda r: r["version"] == "2.7.0" and r["set"] == "both"),
          ("only 2.7.0\nquantified", lambda r: r["version"] == "2.7.0" and r["set"] == "only_2.7.0")]
cols = {"Human": "#4C72B0", "Yeast": "#DD8452", "Ecoli": "#55A868"}
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1.15, 1]})
w = 0.26
for j, s in enumerate(T):
    for i, (lab, f) in enumerate(groups):
        v, n = share([r for r in rows if f(r) and r["species"] == s])
        a1.bar(i + (j - 1) * w, v, w, color=cols[s], alpha=0.85, label=s.replace("Ecoli", "E. coli") if i == 0 else None)
        a1.text(i + (j - 1) * w, v + 1, f"{v:.0f}", ha="center", fontsize=8)
a1.set_xticks(range(len(groups))); a1.set_xticklabels([g[0] for g in groups], fontsize=9); a1.set_ylim(0, 118)
a1.set_ylabel("% within ±0.5 of the true log2 ratio"); a1.legend(frameon=False, loc="upper center", ncol=3)
a1.set_title("Accuracy by protein set (shared = quantified by both versions)", fontsize=10)
e = [r for r in rows if r["version"] == "2.7.0" and r["species"] == "Ecoli"]
for setname, col, lab in [("both", "#9e9e9e", "shared with 2.3.2"), ("only_2.7.0", "#d62728", "quantified only by 2.7.0")]:
    p = [r for r in e if r["set"] == setname]
    a2.scatter([r["intA"] for r in p], [r["fc"] for r in p], s=10, alpha=0.6, color=col, label=f"{lab} (n={len(p)})")
a2.axhline(2, ls="--", color="k", lw=1); a2.axhspan(1.5, 2.5, color="k", alpha=0.06)
a2.set_xlabel("mean log2 intensity in condition A (E. coli at 5%)"); a2.set_ylabel("log2(B/A), truth = +2")
a2.set_title("E. coli proteins, DIA-NN 2.7.0", fontsize=10); a2.legend(frameon=False, loc="lower right")
fig.suptitle("DIA-NN 2.7.0 vs 2.3.2 on PXD028735 (SWATH, 3 vs 3, same settings)", fontsize=11)
fig.tight_layout(); out = sys.argv[1] if len(sys.argv) > 1 else "diann-270-vs-232.png"; fig.savefig(out, dpi=130); print("saved", out)
