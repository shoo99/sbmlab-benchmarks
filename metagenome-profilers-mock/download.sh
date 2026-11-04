#!/usr/bin/env bash
# Downloads the mock-community runs and the Kraken2 databases into data/ and db/ (not tracked).
#   data/ERR2984773.sra   ZymoBIOMICS D6300 (even), Illumina MiSeq 2x151, 8.8 M pairs, 1.5 GB (Nicholls et al. 2019)
#   data/ERR2935805.sra   ZymoBIOMICS D6310 (log distribution), Illumina HiSeq 1500 2x101, 47.8 M pairs, 4.4 GB
#                         (both from NCBI SRA's AWS open-data mirror, s3://sra-pub-run-odp)
#   db/k2_standard_08gb   Kraken2 standard database capped at 8 GB, 2024-09-04 (archaea, bacteria, viral, plasmid,
#                         human, UniVec_Core; no fungi), with Bracken k-mer distributions
#   db/k2_pluspf_08gb     Kraken2 PlusPF capped at 8 GB, 2024-09-04 (standard + protozoa + fungi)
#   (optional) db/k2_standard   the uncapped standard database (about 80 GB on disk and in RAM): DB_FULL=1
# Usage: bash download.sh [D6300|D6310 ...]   (default both; existing files are kept)
set -euo pipefail
cd "$(dirname "$0")" || exit 1
mkdir -p data db
SRA=https://sra-pub-run-odp.s3.amazonaws.com/sra
K2=https://genome-idx.s3.amazonaws.com/kraken
declare -A RUN=([D6300]=ERR2984773 [D6310]=ERR2935805)
SETS=("$@"); [ ${#SETS[@]} -eq 0 ] && SETS=(D6300 D6310)
for s in "${SETS[@]}"; do
  r=${RUN[$s]}
  [ -s "data/$r.sra" ] || { echo "get  $r ($s)"; curl -fsSL --retry 3 -o "data/$r.sra.part" "$SRA/$r/$r" && mv "data/$r.sra.part" "data/$r.sra"; }
done
dbs=(k2_standard_08gb_20240904 k2_pluspf_08gb_20240904); [ "${DB_FULL:-0}" = 1 ] && dbs+=(k2_standard_20240904)
for d in "${dbs[@]}"; do
  name=${d%_2024*}
  [ -s "db/$name/hash.k2d" ] && { echo "have $name"; continue; }
  echo "get  $name (streamed and unpacked)"
  mkdir -p "db/$name" && curl -fsSL --retry 3 "$K2/$d.tar.gz" | tar -xz -C "db/$name"
done
(md5sum data/*.sra; for d in db/*/; do echo "$(md5sum < "$d/hash.k2d" | cut -d' ' -f1)  $d/hash.k2d"; done) | tee data/md5.txt
