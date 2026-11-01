#!/bin/bash
# DIA-NN 2.7.0 on the PXD028735 HYE SWATH set (2026-10-31), same settings as the 2.3.2 run (README in sbmlab-benchmarks/diann-hye).
# Fresh predicted library with 2.7.0, then MBR search. .quant files go to v270/quant so the 2.3.2 .quant files in raw/ stay intact.
cd /mnt/c/diann-hye
D=/mnt/c/DIA-NN/2.7.0/diann.exe
F=$(cat cmd_files.txt)
echo "== lib start $(date)"
/usr/bin/time -v $D --fasta 'C:\diann-hye\hye_sp.fasta' --fasta-search --predictor --gen-spec-lib --out-lib 'C:\diann-hye\v270\lib\hye.parquet' --cut 'K*,R*' --missed-cleavages 1 --min-pep-len 7 --max-pep-len 30 --min-pr-charge 2 --max-pr-charge 4 --min-pr-mz 300 --max-pr-mz 1800 --min-fr-mz 200 --max-fr-mz 1800 --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --threads 24 --verbose 1 > v270/libgen.log 2>&1
echo "== lib end $(date) exit $?"; ls -la v270/lib
L=$(ls v270/lib | grep -E 'predicted\.speclib$' | head -1); [ -z "$L" ] && L=$(ls v270/lib | grep -E '\.parquet$' | head -1)
echo "library file: $L"
echo "== search start $(date)"
/usr/bin/time -v $D $F --lib "C:\\diann-hye\\v270\\lib\\$L" --fasta 'C:\diann-hye\hye_sp.fasta' --out 'C:\diann-hye\v270\out\report.parquet' --temp 'C:\diann-hye\v270\quant' --threads 24 --verbose 1 --qvalue 0.01 --matrices --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --reanalyse --rt-profiling > v270/run.log 2>&1
echo "== search end $(date) exit $?"
