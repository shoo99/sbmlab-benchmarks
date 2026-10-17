#!/bin/bash
# Clean rerun for the quantification-mode comparison (2026-10-16): fixed mass accuracy (MS2 16 ppm, MS1 13 ppm, window 9,
# values from the original run's log), MBR on, --export-quant. Step 1 saves first-pass .quant files; steps 2-3 reuse them.
cd /mnt/c/diann-hye
D=/mnt/c/DIA-NN/2.3.2/diann.exe
F=$(cat cmd_files.txt)
C="--lib C:\\diann-hye\\lib\\hye.predicted.speclib --fasta C:\\diann-hye\\hye_sp.fasta --threads 24 --verbose 1 --qvalue 0.01 --matrices --export-quant --rt-profiling --met-excision --cut K*,R* --missed-cleavages 1 --var-mods 1 --var-mod UniMod:35,15.994915,M --unimod4 --mass-acc 16 --mass-acc-ms1 13 --window 9"
run() { name=$1; shift; mkdir -p c_$name; echo "== $name start $(date)"; /usr/bin/time -v $D $F $C "$@" --out "C:\\diann-hye\\c_$name\\report.parquet" > c_$name/run.log 2>&1; echo "== $name end $(date) exit $?"; }
mkdir -p quant_backup && mv raw/*.quant quant_backup/ 2>/dev/null   # keep the old auto-mass-acc .quant files aside
run mbr_quantums --reanalyse                              # step 1: fresh, MBR, default QuantUMS (writes first-pass .quant)
run nombr_quantums --use-quant                            # step 2: no-MBR variants from the clean .quant
run nombr_highacc --use-quant --high-acc
run nombr_direct --use-quant --direct-quant
run mbr_highacc --use-quant --reanalyse --high-acc        # step 3: MBR variants reusing first-pass .quant (second pass re-reads raw)
run mbr_direct --use-quant --reanalyse --direct-quant
echo "ALL DONE $(date)"
