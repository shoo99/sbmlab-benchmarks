# STRING confidence score and hub proteins, measured (STRING 12.5, human)

Blog post: [STRING confidence score and hub proteins, measured](https://sbmlab.com/blog/string-confidence-score-hub-proteins-measured-2026)

## Inputs (public downloads, 19 October 2026; not copied here)

| File | Source |
|---|---|
| `9606.protein.links.detailed.v12.5.txt.gz`, `9606.protein.info.v12.5.txt.gz` | stringdb-downloads.org (STRING 12.5) |
| `gene2pubmed.gz`, `Homo_sapiens.gene_info.gz` | ftp.ncbi.nlm.nih.gov/gene/DATA |
| `CEGv2.txt`, `NEGv1.txt` | github.com/hart-lab/bagel (Hart et al. 2017; 2014) |
| `rna_celline.tsv.zip` | proteinatlas.org (data file dated 2023-05-24, 1,206 cell lines) |

Put them in one folder and pass it as the first argument to `hubs.py` (default `/home/sysoft/data/string125`).
The detected proteome is `../diann-hye/report.pg_matrix.tsv` from this repository.

## Scripts (run in this order)

| Script | What it does | Output |
|---|---|---|
| `hubs.py [folder]` | Recomputes every combined score from the seven channel scores with STRING's published formula (prior 0.041), and again without text mining; links per cut-off and channel; partners vs PubMed records; core-essential (CEGv2) vs non-essential (NEGv1) and vs all broadly expressed genes (HPA nTPM ≥ 1 in ≥ 90% of cell lines), with expression and publication controls | `hubs_results.json`, `degree_table.tsv.gz`, `hubs.out` |
| `api_lists.py` | 20 random lists of 200 detected proteins and 20 lists of 15 cholesterol-biosynthesis + 185 random detected proteins → STRING API `ppi_enrichment` at 0.4 and 0.7, whole-genome vs detected background; in-list hubs from the local network | `api_results.json`, `api.out` |
| `summarize_api.py` | PubMed records of in-list hubs vs the detected set | `api_summary_extra.json` |
| `figure.py` | Three-panel figure used in the post | `string-hubs-measured.png` (in the blog repository) |

STRING's own reference implementation of the score is `combine_subscores.v2.py` on stringdb-downloads.org.
