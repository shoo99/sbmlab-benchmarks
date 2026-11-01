#!/bin/bash
# Search step only (library already built by run270.sh). First attempt skipped all .wiff: 2.7.0 ships without the SCIEX
# Clearcore/Sciex DLLs; copied the 23 from the 2.3.2 folder (originally from ProteoWizard) on 2026-10-31 21:55.
cd /mnt/c/diann-hye
D=/mnt/c/DIA-NN/2.7.0/diann.exe
F=$(cat cmd_files.txt)
echo "== search start $(date)"
/usr/bin/time -v $D $F --lib 'C:\diann-hye\v270\lib\hye.predicted.speclib' --fasta 'C:\diann-hye\hye_sp.fasta' --out 'C:\diann-hye\v270\out\report.parquet' --temp 'C:\diann-hye\v270\quant' --threads 24 --verbose 1 --qvalue 0.01 --matrices --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --reanalyse --rt-profiling > v270/run.log 2>&1
echo "== search end $(date) exit $?"
