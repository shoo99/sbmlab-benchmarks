#!/usr/bin/env bash
# Builds the search databases for the entrapment searches into db/.
#   target.fasta                 UniProt Swiss-Prot human 9606 + yeast 559292 + E. coli 83333 (the HYE targets).
#                                Better: TARGET_FASTA=/path/to/hye_sp.fasta, the file the earlier HYE searches used
#                                (fetched 2026-10-06), so that only the entrapment differs.
#   arabidopsis.fasta            UniProt Swiss-Prot Arabidopsis thaliana 3702: the foreign entrapment proteome
#   target_plus_entrapment.fasta, entrapment.fasta, entrapment_info.json   from make_entrapment.py
# Usage: bash make_db.sh [foreign|shuffle]    (default foreign; shuffle makes one shuffled copy of each target, r ~ 1)
set -euo pipefail
cd "$(dirname "$0")" || exit 1
MODE=${1:-foreign}
mkdir -p db
UNIPROT="https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28reviewed%3Atrue%29%20AND%20%28organism_id%3A"
if [ -n "${TARGET_FASTA:-}" ]; then
  cp "$TARGET_FASTA" db/target.fasta
elif [ ! -s db/target.fasta ]; then
  : > db/target.fasta.part
  for taxon in 9606 559292 83333; do curl -fsSL --retry 3 "${UNIPROT}${taxon}%29" >> db/target.fasta.part; done
  mv db/target.fasta.part db/target.fasta
fi
if [ "$MODE" = foreign ]; then
  [ -s db/arabidopsis.fasta ] || { curl -fsSL --retry 3 -o db/arabidopsis.fasta.part "${UNIPROT}3702%29" && mv db/arabidopsis.fasta.part db/arabidopsis.fasta; }
  python3 make_entrapment.py db/target.fasta foreign db/arabidopsis.fasta --out db
else
  python3 make_entrapment.py db/target.fasta shuffle 1 --out db
fi
grep -c '^>' db/target.fasta db/entrapment.fasta db/target_plus_entrapment.fasta
(cd db && md5sum ./*.fasta > md5.txt)
echo "next: search the PXD028735 runs against db/target_plus_entrapment.fasta (see README), then python3 extract.py searches.tsv"
