#!/usr/bin/env bash
# Lists the sequencing runs of a study at ENA, to choose the runs for runs.tsv.
# Usage: ./list_runs.sh [study accession] > all_runs.tsv
#   default SRP025982, the SRA study of the SEQC RNA-seq series GEO GSE47792 (check the titles: they name the site,
#   sample A-D, replicate and lane).
set -euo pipefail
STUDY=${1:-SRP025982}
curl -fsS --retry 3 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$STUDY&result=read_run&format=tsv&limit=0&fields=run_accession,sample_title,experiment_title,library_layout,read_count,base_count"
