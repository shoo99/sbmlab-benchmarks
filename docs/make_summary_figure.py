"""Regenerate docs/known-ratio-summary.png and docs/social-preview.png from result tables in this repository.

Run from the repository root:  python docs/make_summary_figure.py
Needs matplotlib (tested with 3.11). Reads only these files:
  maxquant-hye/quant_scores.csv    MaxQuant 2.8.1 LFQ, each engine's own protein list
  fragpipe-hye/fragpipe_scores.csv FragPipe 24 MaxLFQ, own list
  diann-hye/diann_scores.csv       DIA-NN 2.3.2 pg_matrix (MBR, default filter), own list
  diann-hye/common_3866.csv        log2(B/A) of the 3,866 proteins all three engines quantified
  cv-hye/cv_summary.csv            median replicate CV, all species
Scoring rule (same as the folders): a protein counts if log2(B/A) is within 0.5 of the truth
(human 0, yeast -1, E. coli +2). Every printed number is also written to stdout so it can be checked.
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs")
TRUTH = {"Human": 0.0, "Yeast": -1.0, "Ecoli": 2.0}
SPECIES = ["Human", "Yeast", "Ecoli"]
SP_LABEL = {"Human": "Human (1:1)", "Yeast": "Yeast (1:2)", "Ecoli": "E. coli (4:1)"}
ENGINES = ["MaxQuant", "FragPipe", "DIA-NN"]
ENG_LABEL = {"MaxQuant": "MaxQuant 2.8.1 LFQ", "FragPipe": "FragPipe 24 MaxLFQ", "DIA-NN": "DIA-NN 2.3.2 (SWATH)"}
COLOR = {"MaxQuant": "#2a78d6", "FragPipe": "#eb6834", "DIA-NN": "#1baf7a"}  # colourblind-checked order
MARKER = {"MaxQuant": "o", "FragPipe": "s", "DIA-NN": "D"}                    # shape as a second cue
GREY = "#8a8985"
INK, INK2, GRID, BG = "#0b0b0b", "#52514e", "#e4e3df", "#ffffff"


def rows(path):
    with open(os.path.join(ROOT, path), newline="") as f:
        return list(csv.DictReader(f))


def own_lists():
    """Percent within +/-0.5 and n, each engine scored on every protein it quantified."""
    out = {}
    for r in rows("maxquant-hye/quant_scores.csv"):
        if r["quant"] == "LFQ intensity":
            out[("MaxQuant", r["species"])] = (100 * float(r["within_0.5"]), int(r["n"]))
    for r in rows("fragpipe-hye/fragpipe_scores.csv"):
        if r["quant"] == "MaxLFQ Intensity":
            out[("FragPipe", r["species"])] = (100 * float(r["within_0.5"]), int(r["n"]))
    for r in rows("diann-hye/diann_scores.csv"):
        if r["set"].startswith("pg_matrix (MBR"):
            out[("DIA-NN", r["species"])] = (100 * float(r["within_0.5"]), int(r["n"]))
    return out


def common_set():
    """Percent within +/-0.5 on the proteins all three engines quantified (recomputed from log2 ratios)."""
    data = rows("diann-hye/common_3866.csv")
    col = {"MaxQuant": "fc_mq", "FragPipe": "fc_fp", "DIA-NN": "fc_dn"}
    out, n_total = {}, len(data)
    for sp in SPECIES:
        sub = [r for r in data if r["species"] == sp]
        for e in ENGINES:
            hit = sum(abs(float(r[col[e]]) - TRUTH[sp]) <= 0.5 for r in sub)
            out[(e, sp)] = (100 * hit / len(sub), len(sub))
    return out, n_total


def median_cv():
    want = {"MaxQuant LFQ": "MaxQuant", "FragPipe MaxLFQ": "FragPipe", "DIA-NN pg_matrix": "DIA-NN",
            "MaxQuant raw Intensity": "MaxQuant raw"}
    return {want[r["engine"]]: float(r["median_cv"]) for r in rows("cv-hye/cv_summary.csv")
            if r["species"] == "All" and r["engine"] in want}


def style(ax):
    ax.set_facecolor(BG)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(INK2)
    ax.tick_params(colors=INK2, length=0, labelsize=12)
    ax.grid(axis="x", color=GRID, lw=1)
    ax.set_axisbelow(True)


def dot_panel(ax, vals, title, show_species=True):
    style(ax)
    offs = {"MaxQuant": 0.24, "FragPipe": 0.0, "DIA-NN": -0.24}
    for i, sp in enumerate(SPECIES):
        y0 = len(SPECIES) - 1 - i
        for e in ENGINES:
            v, _ = vals[(e, sp)]
            y = y0 + offs[e]
            ax.plot([70, v], [y, y], color=GRID, lw=1.2, zorder=1, solid_capstyle="round")
            ax.scatter([v], [y], s=70, marker=MARKER[e], color=COLOR[e], edgecolor=BG, linewidth=1.5, zorder=3)
            ax.text(v + 0.9, y, f"{v:.1f}", va="center", ha="left", fontsize=11.5, color=INK)
    ax.set_xlim(70, 104.5)
    ax.set_xticks([70, 80, 90, 100])
    ax.set_xticklabels(["70%", "80%", "90%", "100%"])
    ax.set_ylim(-0.6, 2.6)
    ax.set_yticks([2, 1, 0])
    ax.set_yticklabels([SP_LABEL[s] for s in SPECIES] if show_species else [], fontsize=12.5, color=INK)
    ax.set_title(title, loc="left", fontsize=13, color=INK, pad=10)


def summary_figure(own, common, n_common, cv):
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig = plt.figure(figsize=(13, 5.6), dpi=170, facecolor=BG)
    gs = fig.add_gridspec(1, 2, wspace=0.10, left=0.105, right=0.665, top=0.80, bottom=0.25)
    a, b = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    c = fig.add_subplot(fig.add_gridspec(1, 1, left=0.79, right=0.985, top=0.80, bottom=0.25)[0])
    dot_panel(a, own, "A  Each engine's own protein list")
    dot_panel(b, common, f"B  The {n_common:,} proteins all three quantified", show_species=False)

    style(c)
    order = ["DIA-NN", "MaxQuant", "FragPipe", "MaxQuant raw"]
    labels = {"DIA-NN": "DIA-NN", "MaxQuant": "MaxQuant LFQ", "FragPipe": "FragPipe", "MaxQuant raw": "MaxQuant raw\nIntensity"}
    for i, k in enumerate(order):
        y = len(order) - 1 - i
        col = COLOR.get(k, GREY)
        c.barh(y, cv[k], height=0.56, color=col, edgecolor=BG, linewidth=2)
        c.text(cv[k] + 0.4, y, f"{cv[k]:.1f}%", va="center", fontsize=11.5, color=INK)
    c.set_yticks(range(len(order) - 1, -1, -1))
    c.set_yticklabels([labels[k] for k in order], fontsize=12, color=INK)
    c.set_xlim(0, 26)
    c.set_xticks([0, 10, 20])
    c.set_xticklabels(["0%", "10%", "20%"])
    c.set_title("C  Median replicate CV", loc="left", fontsize=13, color=INK, pad=10)

    fig.text(0.012, 0.955, "Proteins within ±0.5 log2 of the true B/A ratio, and replicate CV",
             fontsize=16, color=INK, weight="bold", ha="left")
    fig.text(0.012, 0.895, "PRIDE PXD028735. True B:A ratio by design: human 1:1, yeast 1:2, E. coli 4:1. Three runs per condition, match-between-runs on.",
             fontsize=12, color=INK2, ha="left")
    handles = [plt.Line2D([], [], ls="", marker=MARKER[e], markersize=9, color=COLOR[e], label=ENG_LABEL[e]) for e in ENGINES]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.10, 0.115), ncol=3, frameon=False, fontsize=12,
               handletextpad=0.4, columnspacing=1.6)
    fig.text(0.012, 0.03,
             "MaxQuant and FragPipe: the same 6 Orbitrap DDA runs. DIA-NN: SWATH runs of the same samples on a TripleTOF 6600+ (different instrument).\n"
             "Sources: A maxquant-hye/quant_scores.csv, fragpipe-hye/fragpipe_scores.csv, diann-hye/diann_scores.csv; "
             "B diann-hye/common_3866.csv; C cv-hye/cv_summary.csv.",
             fontsize=10, color=INK2, ha="left", linespacing=1.45)
    path = os.path.join(OUT, "known-ratio-summary.png")
    fig.savefig(path, facecolor=BG)
    plt.close(fig)
    return path


def social_preview(own, common, cv):
    """1280 x 640 image for GitHub's social preview (Settings > Social preview)."""
    fig = plt.figure(figsize=(12.8, 6.4), dpi=100, facecolor=BG)
    fig.text(0.06, 0.82, "sbmlab-benchmarks", fontsize=46, weight="bold", color=INK)
    fig.text(0.06, 0.72, "MaxQuant, FragPipe and DIA-NN scored against known ratios", fontsize=24, color=INK2)
    fig.text(0.06, 0.615, "E. coli proteins within ±0.5 log2 of the true 4:1 ratio, on the "
             f"{common[('MaxQuant', 'Ecoli')][1]} E. coli proteins all three quantified", fontsize=15, color=INK2)
    x = 0.06
    for e in ENGINES:
        fig.text(x, 0.43, f"{common[(e, 'Ecoli')][0]:.1f}%", fontsize=52, weight="bold", color=INK)
        fig.patches.append(matplotlib.patches.Rectangle((x, 0.395), 0.035, 0.012, transform=fig.transFigure, color=COLOR[e]))
        fig.text(x, 0.33, ENG_LABEL[e], fontsize=17, color=INK)
        fig.text(x, 0.28, f"median replicate CV {cv[e]:.1f}%", fontsize=15, color=INK2)
        x += 0.31
    fig.text(0.06, 0.13, "Public HYE benchmark PRIDE PXD028735 · scripts, settings and result tables · MIT · DOI 10.5281/zenodo.23355028",
             fontsize=14, color=INK2)
    fig.text(0.06, 0.06, "sbmlab.com", fontsize=20, weight="bold", color=COLOR["MaxQuant"])
    path = os.path.join(OUT, "social-preview.png")
    fig.savefig(path, facecolor=BG)
    plt.close(fig)
    return path


if __name__ == "__main__":
    own, (common, n_common), cv = own_lists(), common_set(), median_cv()
    for name, vals in (("own lists", own), (f"common {n_common}", common)):
        for e in ENGINES:
            print(f"{name:12s} {e:9s} " + "  ".join(f"{s} {vals[(e, s)][0]:.1f}% (n {vals[(e, s)][1]})" for s in SPECIES))
    print("median CV  " + "  ".join(f"{k} {v:.1f}%" for k, v in cv.items()))
    print(summary_figure(own, common, n_common, cv))
    print(social_preview(own, common, cv))
