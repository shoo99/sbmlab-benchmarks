"""Scores the GO release drift study from the tables enrich.py wrote to results/ (issue #19).

Per release and evidence setting (all codes; without IEA):
  random lists (20)   false_terms: significant terms summed over the 20 lists (truth: 0)
                      lists_with_false_terms
  planted lists (5)   planted_recovery: share of lists with "cholesterol biosynthetic process" (GO:0006695) at
                      FDR < 0.05; planted_median_rank: its median rank by p-value among the tested terms;
                      terms_found: mean significant terms per list
  cholesterol genes   genes_annotated: of the detected Reactome cholesterol-biosynthesis proteins, how many carry
                      GO:0006695 (directly or through a descendant)
  airway up / down    terms_found
  airway up           fdr and rank of "response to glucocorticoid" (GO:0051384) and "cellular response to
                      glucocorticoid stimulus" (GO:0071385): dexamethasone is a glucocorticoid, so these are the
                      expected terms of the real list
  every list group    jaccard_vs_2025, jaccard_vs_previous: Jaccard overlap of the significant terms with those of
                      the 2025-10-10 release / of the previous release (terms renamed or merged are compared by
                      their 2025 id; terms obsolete by 2025 count as lost); terms_obsolete_by_2025
  backgrounds         annotated (members with at least one term), mean_terms_per_member, terms_tested
Ontology vs annotations: one run reads the 2016 annotations through the 2025 ontology (gaf2016+obo2025);
  jaccard_vs_2016 then measures the change of ontology alone and jaccard_vs_2025 the change of annotations alone.
Between releases (detected proteome): proteins gained or lost (annotation present in one release only),
share of proteins with an unchanged propagated term set, mean per-protein Jaccard of the propagated term sets,
mean per-term Jaccard of the member sets of terms tested in both releases.
Usage: python3 score.py [results_dir] [out_dir]   (standard library only)
"""
import csv
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
res = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "results"
out = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE
FOLDER = "go-release-drift"
COLUMNS = ["folder", "dataset", "tool", "version", "setting", "subset", "metric", "value", "n"]
DATASET = "unchanged gene lists (DIA-NN HYE human proteome; airway dexamethasone DE) vs GO releases 2016-2025"
TOOL = "GO over-representation (hypergeometric, BH)"
RELEASES = ["2016-01-01", "2020-01-01", "2023-01-01", "2025-10-10"]
HYBRID, HYBRID_LABEL = "gaf2016+obo2025", "GO 2016-01-01 annotations, 2025-10-10 ontology"
EVID = {"all": "all evidence codes", "no_IEA": "without IEA"}
FDR = 0.05


def read(name):
    return list(csv.DictReader(open(res / name), delimiter="\t"))


terms, lists, tracked = read("terms.tsv"), read("lists.tsv"), read("tracked.tsv")
planted = [r for r in tracked if r["list"].startswith("planted")]
GLUCO = {"GO:0051384": "response to glucocorticoid", "GO:0071385": "cellular response to glucocorticoid stimulus"}
annot, changes = read("annotation.tsv"), read("changes.tsv")
groups = {"random lists (20)": [f"random_{i}" for i in range(1, 21)],
          "planted lists (5)": [f"planted_{i}" for i in range(1, 6)],
          "airway up": ["airway_up"], "airway down": ["airway_down"]}
sig = {}  # (list, release, evidence) -> set of significant terms by their 2025 id
for r in terms:
    key = (r["list"], r["release"], r["evidence"])
    sig.setdefault(key, set()).add(r["term_in_2025"] if r["term_in_2025"] != "obsolete" else "obsolete:" + r["term"])
nsig = {(r["list"], r["release"], r["evidence"]): int(r["terms_significant"]) for r in lists}


def jac(a, b):
    return len(a & b) / len(a | b) if a | b else 1.0


rows = []


def add(release, ev, subset, metric, value, n):
    rows.append([FOLDER, DATASET, TOOL, f"GO {release}", EVID.get(ev, ev), subset, metric, value, n])


print(f"== {DATASET}\n   {TOOL}; FDR < {FDR}; terms with 10-500 annotated background members\n")
for ev in EVID:
    print(f"-- {EVID[ev]}")
    print(f"   {'release':<11} {'random: false terms':>20} {'lists':>6} {'planted found':>14} {'median rank':>12} "
          f"{'chol. genes':>12} {'airway up':>10} {'airway down':>12} {'glucocorticoid FDR (rank)':>27}")
    for i, rel in enumerate(RELEASES):
        rnd = groups["random lists (20)"]
        false_terms = sum(nsig[(l, rel, ev)] for l in rnd)
        lists_false = sum(nsig[(l, rel, ev)] > 0 for l in rnd)
        add(rel, ev, "random lists (20)", "false_terms", false_terms, 20)
        add(rel, ev, "random lists (20)", "lists_with_false_terms", lists_false, 20)
        pl = [r for r in planted if r["release"] == rel and r["evidence"] == ev and r["term"] == "GO:0006695"]
        found = [r["fdr"] != "" and float(r["fdr"]) < FDR for r in pl]
        ranks = [int(r["rank"]) for r in pl if r["rank"]]
        recovery = sum(found) / len(pl)
        med_rank = statistics.median(ranks) if ranks else ""
        add(rel, ev, "planted lists (5)", "planted_recovery", round(recovery, 4), len(pl))
        add(rel, ev, "planted lists (5)", "planted_median_rank", med_rank, len(ranks))
        found_terms = statistics.mean(nsig[(l, rel, ev)] for l in groups["planted lists (5)"])
        add(rel, ev, "planted lists (5)", "terms_found", round(found_terms, 1), 5)
        a = next(r for r in annot if r["release"] == rel and r["evidence"] == ev and r["background"] == "proteins")
        chol = int(a["cholesterol_genes_with_GO0006695"])
        add(rel, ev, "detected cholesterol-biosynthesis proteins (16)", "genes_annotated", chol, 16)
        for g in ("airway up", "airway down"):
            l = groups[g][0]
            add(rel, ev, g, "terms_found", nsig[(l, rel, ev)], 1)
        gl = {r["term"]: r for r in tracked if r["list"] == "airway_up" and r["release"] == rel and r["evidence"] == ev}
        for t, name in GLUCO.items():
            if gl[t]["fdr"]:                       # not tested when fewer than 10 background genes carry the term
                add(rel, ev, "airway up", f"fdr: {name}", float(gl[t]["fdr"]), 1)
                add(rel, ev, "airway up", f"rank: {name}", int(gl[t]["rank"]), 1)
        gtxt = f"{float(gl['GO:0051384']['fdr']):.3f} ({gl['GO:0051384']['rank']})" if gl["GO:0051384"]["fdr"] else "not tested"
        for g, members in groups.items():
            if g.startswith("random"):
                continue
            j_new = statistics.mean(jac(sig.get((l, rel, ev), set()), sig.get((l, RELEASES[-1], ev), set()))
                                    for l in members)
            obsolete = statistics.mean(sum(t.startswith("obsolete:") for t in sig.get((l, rel, ev), set()))
                                       for l in members)
            if rel != RELEASES[-1]:
                add(rel, ev, g, "jaccard_vs_2025", round(j_new, 4), len(members))
                add(rel, ev, g, "terms_obsolete_by_2025", round(obsolete, 1), len(members))
            if i > 0:
                j_prev = statistics.mean(jac(sig.get((l, rel, ev), set()), sig.get((l, RELEASES[i - 1], ev), set()))
                                         for l in members)
                add(rel, ev, g, "jaccard_vs_previous", round(j_prev, 4), len(members))
        for bg, label in (("proteins", "detected proteome (5,661 proteins)"), ("genes", "airway tested genes (13,195)")):
            b = next(r for r in annot if r["release"] == rel and r["evidence"] == ev and r["background"] == bg)
            add(rel, ev, label, "annotated", int(b["annotated"]), int(b["units"]))
            add(rel, ev, label, "mean_terms_per_member", float(b["mean_terms_per_unit"]), int(b["annotated"]))
            add(rel, ev, label, "terms_tested", int(b["terms_tested"]), int(b["annotated"]))
        print(f"   {rel:<11} {false_terms:>20} {lists_false:>6} {sum(found):>9} of {len(pl)} {med_rank:>12} "
              f"{chol:>7} of 16 {nsig[('airway_up', rel, ev)]:>10} {nsig[('airway_down', rel, ev)]:>12} {gtxt:>27}")
    print(f"\n   Jaccard of the significant terms with 2025-10-10 (planted: mean of 5 lists; obsolete terms by 2025)")
    print(f"   {'release':<11} {'planted':>8} {'airway up':>10} {'airway down':>12} {'obsolete, up/down':>18}")
    for rel in RELEASES[:-1]:
        js = [statistics.mean(jac(sig.get((l, rel, ev), set()), sig.get((l, RELEASES[-1], ev), set())) for l in m)
              for g, m in groups.items() if not g.startswith("random")]
        ob = [sum(t.startswith("obsolete:") for t in sig.get((l, rel, ev), set())) for l in ("airway_up", "airway_down")]
        print(f"   {rel:<11} {js[0]:>8.3f} {js[1]:>10.3f} {js[2]:>12.3f} {ob[0]:>11} / {ob[1]}")
    print(f"\n   ontology vs annotations: 2016 annotations read through the 2025 ontology (hybrid)")
    print(f"   {'list':<18} {'terms 2016':>10} {'hybrid':>7} {'2025':>5} {'Jaccard 2016-hybrid':>20} {'hybrid-2025':>12} "
          f"{'2016-2025':>10}")
    for g, members in groups.items():
        if g.startswith("random"):
            continue
        j = {}
        for a, b in ((RELEASES[0], HYBRID), (HYBRID, RELEASES[-1]), (RELEASES[0], RELEASES[-1])):
            j[(a, b)] = statistics.mean(jac(sig.get((l, a, ev), set()), sig.get((l, b, ev), set())) for l in members)
        nt = {r: statistics.mean(nsig[(l, r, ev)] for l in members) for r in (RELEASES[0], HYBRID, RELEASES[-1])}
        rows.append([FOLDER, DATASET, TOOL, HYBRID_LABEL, EVID[ev], g, "terms_found", round(nt[HYBRID], 1), len(members)])
        rows.append([FOLDER, DATASET, TOOL, HYBRID_LABEL, EVID[ev], g, "jaccard_vs_2016",
                     round(j[(RELEASES[0], HYBRID)], 4), len(members)])
        rows.append([FOLDER, DATASET, TOOL, HYBRID_LABEL, EVID[ev], g, "jaccard_vs_2025",
                     round(j[(HYBRID, RELEASES[-1])], 4), len(members)])
        print(f"   {g:<18} {nt[RELEASES[0]]:>10.1f} {nt[HYBRID]:>7.1f} {nt[RELEASES[-1]]:>5.1f} "
              f"{j[(RELEASES[0], HYBRID)]:>20.3f} {j[(HYBRID, RELEASES[-1])]:>12.3f} {j[(RELEASES[0], RELEASES[-1])]:>10.3f}")
    rnd_h = sum(nsig[(l, HYBRID, ev)] for l in groups["random lists (20)"])
    rows.append([FOLDER, DATASET, TOOL, HYBRID_LABEL, EVID[ev], "random lists (20)", "false_terms", rnd_h, 20])
    print(f"   random lists: {rnd_h} false terms with the hybrid")
    print(f"\n   between releases, detected proteome: annotation of each protein and members of each tested term")
    print(f"   {'releases':<30} {'gained':>7} {'lost':>5} {'unchanged':>10} {'protein Jaccard':>16} {'term Jaccard':>13}")
    for c in changes:
        if c["evidence"] != ev:
            continue
        pair = f"GO {c['release_a']} -> {c['release_b']}".replace(HYBRID, "2016 annotations with 2025 ontology")
        both = int(c["proteins_in_both"])
        unchanged = int(c["proteins_unchanged"]) / both
        rows.append([FOLDER, DATASET, TOOL, pair, EVID[ev], "detected proteome (5,661 proteins)", "proteins_gained",
                     int(c["proteins_gained"]), 5661])
        rows.append([FOLDER, DATASET, TOOL, pair, EVID[ev], "detected proteome (5,661 proteins)", "proteins_lost",
                     int(c["proteins_lost"]), 5661])
        rows.append([FOLDER, DATASET, TOOL, pair, EVID[ev], "detected proteome (5,661 proteins)", "unchanged_share",
                     round(unchanged, 4), both])
        rows.append([FOLDER, DATASET, TOOL, pair, EVID[ev], "detected proteome (5,661 proteins)",
                     "mean_protein_jaccard", float(c["mean_protein_jaccard"]), both])
        rows.append([FOLDER, DATASET, TOOL, pair, EVID[ev], "tested GO terms", "mean_term_jaccard",
                     float(c["mean_term_jaccard"]), int(c["terms_tested_both"])])
        names = f"{c['release_a']} -> {c['release_b']}".replace(HYBRID, "hybrid")
        print(f"   {names:<30} {c['proteins_gained']:>7} {c['proteins_lost']:>5} "
              f"{unchanged:>10.3f} {float(c['mean_protein_jaccard']):>16.3f} {float(c['mean_term_jaccard']):>13.3f}")
    print()
out.mkdir(parents=True, exist_ok=True)
with open(out / "scores.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(COLUMNS)
    w.writerows(rows)
print(f"wrote scores.csv ({len(rows)} rows)")
