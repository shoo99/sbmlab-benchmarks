"""Pig -> mouse homologs from Ensembl BioMart via the plain martservice XML API (what pybiomart/the web UI use).
Run for release 111 (Jan 2024 archive, the tutorial's pin) and 116 (current, Oct 2026).
www.ensembl.org/biomart now answers with HTTP 308 to jun2026.archive.ensembl.org, which biomaRt 2.58 does not follow."""
import urllib.request, urllib.parse, csv, io, collections, sys
HOSTS = {"111": "https://jan2024.archive.ensembl.org", "116": "https://jun2026.archive.ensembl.org"}
ATTR = ["ensembl_gene_id", "external_gene_name", "mmusculus_homolog_ensembl_gene",
        "mmusculus_homolog_associated_gene_name", "mmusculus_homolog_orthology_type"]
xml = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1">'
       '<Dataset name="sscrofa_gene_ensembl" interface="default">' + "".join(f'<Attribute name="{a}"/>' for a in ATTR) + '</Dataset></Query>')
res = {}
for rel, host in HOSTS.items():
    url = host + "/biomart/martservice?query=" + urllib.parse.quote(xml)
    txt = urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "sbmlab-check"}), timeout=600).read().decode()
    rows = list(csv.reader(io.StringIO(txt), delimiter="\t"))
    hdr, rows = rows[0], [r for r in rows[1:] if len(r) == len(ATTR) and r[2]]
    types = collections.Counter(r[4] for r in rows)
    one = [r for r in rows if r[4] == "ortholog_one2one"]
    res[rel] = {(r[0], r[2]) for r in one}
    print(f"Ensembl {rel}: header {hdr[:3]}... rows with mouse homolog {len(rows):,}; types {dict(types)}")
    print(f"  one2one pairs {len(one):,}; pig genes with a name {sum(1 for r in one if r[1]):,}")
    with open(f"one2one_r{rel}.tsv", "w") as f:
        f.write("\t".join(ATTR) + "\n"); f.writelines("\t".join(r) + "\n" for r in one)
print(f"one2one pairs in both releases: {len(res['111'] & res['116']):,}; only 111 {len(res['111'] - res['116']):,}; only 116 {len(res['116'] - res['111']):,}")
