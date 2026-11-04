"""The same gene lists, enriched against GO releases years apart (issue #19).

Releases: 2016-01-01, 2020-01-01, 2023-01-01, 2025-10-10 (go-basic.obo and goa_human.gaf.gz of each, from download.sh).
Gene lists, all fixed before any release is read:
  random_1..20    300 proteins drawn from the 5,661 human protein groups that DIA-NN detected in ../diann-hye
                  (seeds 1-20); no term should be enriched against the detected background
  planted_1..5    15 detected proteins of Reactome "Cholesterol biosynthesis" (R-HSA-191273, as in MSigDB 7.5.1;
                  reactome_cholesterol.tsv) plus 285 other detected proteins (seeds 101-105); the GO term
                  "cholesterol biosynthetic process" (GO:0006695) should be enriched
  airway_up, airway_down   genes up or down with dexamethasone in airway smooth-muscle cells (airway_de.R: DESeq2,
                  padj < 0.05, |log2FC| >= 1); background: all 16,596 tested genes
Identifiers are held fixed: proteins by UniProt accession (the leading accession of each protein group), genes by
Ensembl id, mapped to UniProt accessions with one HGNC table (2025-10-07). A gene carries a term if any of its
accessions does.
Annotation, per release: GAF rows of UniProtKB human entries, without NOT qualifiers; terms propagated to their
ancestors over is_a and part_of in that release's go-basic.obo; alternative ids mapped to the primary id; two
evidence settings, all codes and all codes but IEA (electronic annotations).
Test: one-sided hypergeometric (over-representation) for every term with 10-500 annotated background members, all
three GO aspects together, Benjamini-Hochberg over the tested terms; significant at FDR < 0.05. The background is
the list's universe (detected proteins or tested genes) restricted to members with at least one annotation in that
release, as most ORA tools do.
One extra run reads the 2016 annotations through the 2025 ontology (alternative ids mapped, terms obsolete by 2025
dropped), so that 2016 -> hybrid is the change of ontology alone and hybrid -> 2025 the change of annotations alone.
Writes results/: terms.tsv (significant terms), lists.tsv (list sizes and term counts), tracked.tsv (the sterol
terms in the planted lists, the glucocorticoid-response terms in the airway lists), annotation.tsv (background annotation per release), changes.tsv (per-protein and
per-term changes between releases), cholesterol_genes.tsv (the sterol terms of each planted protein).
Standard library only; under a minute.
Usage: python3 enrich.py [data_dir]
"""
import csv
import gzip
import random
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from go_ora import ancestors_fn, bh, hyper_upper, log_factorials, read_gaf, read_obo  # noqa: E402
data = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "data"
res = HERE / "results"
RELEASES = ["2016-01-01", "2020-01-01", "2023-01-01", "2025-10-10"]
HYBRID = "gaf2016+obo2025"     # 2016 annotations read through the 2025 ontology: separates ontology from annotation change
RUNS = [(r, r, r) for r in RELEASES] + [(HYBRID, "2025-10-10", "2016-01-01")]    # (label, ontology, annotations)
EVIDENCE = {"all": lambda ev: True, "no_IEA": lambda ev: ev != "IEA"}
MIN_SIZE, MAX_SIZE, FDR = 10, 500, 0.05
TRACKED = {"planted": ["GO:0006695", "GO:0016126"],    # cholesterol / sterol biosynthetic process (planted answer)
           "airway": ["GO:0051384", "GO:0071385"]}     # (cellular) response to glucocorticoid: dexamethasone is one
t0 = time.time()


# ---- fixed inputs: detected proteins, Reactome genes, airway genes, HGNC mapping, lists
det, det_gene = [], {}
with open(HERE.parent / "diann-hye" / "report.pg_matrix.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        names = r["Protein.Names"]
        if "_HUMAN" in names and not re.search("_YEAST|_ECOLI", names):
            acc = r["Protein.Group"].split(";")[0]
            if acc not in det_gene:
                det.append(acc)
                det_gene[acc] = r["Genes"].split(";")[0]
det = sorted(det)
chol_sym = [r["symbol"] for r in csv.DictReader(open(HERE / "reactome_cholesterol.tsv"), delimiter="\t")]
chol_det = sorted(a for a in det if det_gene[a] in chol_sym)
hgnc = {}
with open(data / "hgnc_complete_set_2025-10-07.txt", encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["ensembl_gene_id"] and r["uniprot_ids"]:
            hgnc[r["ensembl_gene_id"]] = (r["symbol"], [u.strip() for u in r["uniprot_ids"].split("|") if u.strip()])
with gzip.open(HERE / "airway_de.tsv.gz", "rt") as fh:
    de = [r for r in csv.DictReader(fh, delimiter="\t")]
genes = sorted(r["ensembl"] for r in de if r["ensembl"] in hgnc)
up = {r["ensembl"] for r in de if r["padj"] != "NA" and float(r["padj"]) < 0.05 and float(r["log2fc"]) >= 1}
down = {r["ensembl"] for r in de if r["padj"] != "NA" and float(r["padj"]) < 0.05 and float(r["log2fc"]) <= -1}
print(f"== detected human proteins {len(det)}; Reactome cholesterol-biosynthesis genes {len(chol_sym)}, detected "
      f"{len(chol_det)}; airway genes tested {len(de)}, with UniProt accessions {len(genes)} "
      f"(up {len(up & set(genes))} of {len(up)}, down {len(down & set(genes))} of {len(down)})", flush=True)

lists = {}  # name -> (universe name, members)
for seed in range(1, 21):
    lists[f"random_{seed}"] = ("proteins", random.Random(seed).sample(det, 300))
rest = [a for a in det if a not in set(chol_det)]
for seed in range(1, 6):
    rng = random.Random(100 + seed)
    lists[f"planted_{seed}"] = ("proteins", rng.sample(chol_det, 15) + rng.sample(rest, 285))
lists["airway_up"] = ("genes", sorted(up & set(genes)))
lists["airway_down"] = ("genes", sorted(down & set(genes)))
universes = {"proteins": {a: [a] for a in det}, "genes": {g: hgnc[g][1] for g in genes}}

# ---- per release
terms_rows, list_rows, tracked_rows, annot_rows, chol_rows = [], [], [], [], []
newest_terms, newest_alt = read_obo(data / RELEASES[-1] / "go-basic.obo")
prop_by = {}      # (release, evidence) -> {protein accession: frozenset of terms} for the detected proteins
tested_by = {}    # (release, evidence) -> {term: frozenset of detected proteins}, tested terms only
lf = log_factorials(20000)
for rel, obo_rel, gaf_rel in RUNS:
    terms, alt = (newest_terms, newest_alt) if obo_rel == RELEASES[-1] else read_obo(data / obo_rel / "go-basic.obo")
    anc = ancestors_fn(terms)
    direct = read_gaf(data / gaf_rel / "goa_human.gaf.gz", terms, alt, EVIDENCE)
    print(f"-- {rel}: {sum(not t['obsolete'] for t in terms.values())} terms; annotated accessions "
          f"{len(direct['all'])} (all evidence), {len(direct['no_IEA'])} (no IEA); {time.time() - t0:.0f} s", flush=True)
    for ev in EVIDENCE:
        dir_ev = direct[ev]
        prop_acc = {}

        def prop(acc):
            if acc not in prop_acc:
                s = set()
                for go in dir_ev.get(acc, ()):
                    s |= anc(go)
                prop_acc[acc] = frozenset(s)
            return prop_acc[acc]
        for uname, umap in universes.items():
            unit_terms = {}
            for u, accs in umap.items():
                s = set()
                for a in accs:
                    s |= prop(a)
                if s:
                    unit_terms[u] = s
            N = len(unit_terms)
            members = {}
            for u, s in unit_terms.items():
                for t in s:
                    members.setdefault(t, []).append(u)
            tested = sorted(t for t, m in members.items() if MIN_SIZE <= len(m) <= MAX_SIZE)
            tested_set = set(tested)
            n_direct = sum(len(dir_ev.get(a, ())) for u in unit_terms for a in umap[u])
            annot_rows.append([rel, ev, uname, len(umap), N, n_direct,
                               round(sum(len(s) for s in unit_terms.values()) / max(N, 1), 1), len(tested)])
            if uname == "proteins":
                prop_by[(rel, ev)] = {u: frozenset(s) for u, s in unit_terms.items()}
                tested_by[(rel, ev)] = {t: frozenset(members[t]) for t in tested}
                chol_ann = sum("GO:0006695" in unit_terms.get(a, ()) for a in chol_det)
                for a in chol_det:
                    s = unit_terms.get(a, set())
                    chol_rows.append([det_gene[a], a, rel, ev, int("GO:0006695" in s), int("GO:0016126" in s)])
                annot_rows[-1].append(chol_ann)
            else:
                annot_rows[-1].append("")
            for lname, (lu, lmem) in lists.items():
                if lu != uname:
                    continue
                L = [u for u in lmem if u in unit_terms]
                n = len(L)
                counts = {}
                for u in L:
                    for t in unit_terms[u]:
                        if t in tested_set:
                            counts[t] = counts.get(t, 0) + 1
                p = [hyper_upper(counts.get(t, 0), n, len(members[t]), N, lf) for t in tested]
                q = bh(p)
                order = sorted(range(len(tested)), key=lambda i: (p[i], tested[i]))
                rank = {tested[i]: r + 1 for r, i in enumerate(order)}
                nsig = 0
                for i in order:
                    if q[i] >= FDR:
                        break
                    t = tested[i]
                    nsig += 1
                    nt = newest_alt.get(t, t)
                    nt = nt if nt in newest_terms and not newest_terms[nt]["obsolete"] else "obsolete"
                    terms_rows.append([lname, rel, ev, t, terms[t]["ns"], terms[t]["name"], counts.get(t, 0), n,
                                       len(members[t]), N, f"{p[i]:.4g}", f"{q[i]:.4g}", nt])
                list_rows.append([lname, rel, ev, len(lmem), n, N, len(tested), nsig])
                for t in TRACKED.get(lname.split("_")[0], []):
                    if t in tested_set:
                        i = tested.index(t)
                        tracked_rows.append([lname, rel, ev, t, counts.get(t, 0), len(members[t]), n, N,
                                             f"{p[i]:.4g}", f"{q[i]:.4g}", rank[t], len(tested)])
                    else:
                        tracked_rows.append([lname, rel, ev, t, "", len(members.get(t, ())), n, N, "", "", "",
                                             len(tested)])
    print(f"   enriched {len(lists)} lists x {len(EVIDENCE)} evidence settings; {time.time() - t0:.0f} s", flush=True)

# ---- changes between releases, detected proteins
change_rows = []
pairs = [(RELEASES[i], RELEASES[i + 1]) for i in range(len(RELEASES) - 1)] + [(RELEASES[0], RELEASES[-1]),
                                                                              (RELEASES[0], HYBRID), (HYBRID, RELEASES[-1])]
for ev in EVIDENCE:
    for a, b in pairs:
        pa, pb = prop_by[(a, ev)], prop_by[(b, ev)]
        both = set(pa) & set(pb)
        jac = [len(pa[u] & pb[u]) / len(pa[u] | pb[u]) for u in both]
        same = sum(pa[u] == pb[u] for u in both)
        gained_terms = sum(len(pb[u] - pa[u]) for u in both)
        lost_terms = sum(len(pa[u] - pb[u]) for u in both)
        ta, tb = tested_by[(a, ev)], tested_by[(b, ev)]
        common = set(ta) & set(tb)
        tj = [len(ta[t] & tb[t]) / len(ta[t] | tb[t]) for t in common]
        change_rows.append([a, b, ev, len(set(pb) - set(pa)), len(set(pa) - set(pb)), len(both), same,
                            round(sum(jac) / max(len(jac), 1), 4), gained_terms, lost_terms,
                            len(ta), len(tb), len(common), sum(ta[t] == tb[t] for t in common),
                            round(sum(tj) / max(len(tj), 1), 4)])

res.mkdir(exist_ok=True)


def write(name, header, rows):
    with open(res / name, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


write("terms.tsv", ["list", "release", "evidence", "term", "namespace", "name", "k", "n", "K", "N", "p", "fdr",
                    "term_in_2025"], terms_rows)
write("lists.tsv", ["list", "release", "evidence", "list_size", "n_annotated", "background_annotated",
                    "terms_tested", "terms_significant"], list_rows)
write("tracked.tsv", ["list", "release", "evidence", "term", "k", "K", "n", "N", "p", "fdr", "rank",
                      "terms_tested"], tracked_rows)
write("cholesterol_genes.tsv", ["gene", "accession", "release", "evidence", "GO0006695", "GO0016126"], chol_rows)
write("annotation.tsv", ["release", "evidence", "background", "units", "annotated", "direct_annotations",
                         "mean_terms_per_unit", "terms_tested", "cholesterol_genes_with_GO0006695"], annot_rows)
write("changes.tsv", ["release_a", "release_b", "evidence", "proteins_gained", "proteins_lost", "proteins_in_both",
                      "proteins_unchanged", "mean_protein_jaccard", "protein_terms_gained", "protein_terms_lost",
                      "terms_tested_a", "terms_tested_b", "terms_tested_both", "terms_same_members",
                      "mean_term_jaccard"], change_rows)
print(f"wrote results/ ({len(terms_rows)} significant terms); {time.time() - t0:.0f} s")
print(f"Python {sys.version.split()[0]}")
