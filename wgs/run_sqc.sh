#!/bin/bash
# Scaffold + QC both trizona haplotypes (asm25) and QC the deposited mSicTri.1 for comparison.
# meryl db of all HiFi reads first (shared by the three Merqury runs).
set -euo pipefail
G=/home/toki/gof
O=$G/sqc
mkdir -p $O
cd $O || exit 1
[ "$(pwd)" = "$O" ] || exit 1
export PATH=/home/toki/miniforge3/envs/asm25/bin:$PATH
export TMPDIR=/home/toki/tmp
export MERQURY=$(ls -d /home/toki/miniforge3/envs/asm25/share/merqury*)
step() { echo "#### $* -- $(date '+%F %T')"; }
L=$G/compleasm_lib
H=$G/hic
step "meryl k=21 on HiFi reads"
[ -d hifi.meryl ] || meryl count k=21 threads=64 memory=200 $G/reads/ERR16751478.fastq.gz output hifi.meryl > meryl.log 2>&1
step "hap1 + hap2 scaffold/QC (parallel)"
bash $G/code/wgs/scaffold_qc.sh $G/asm25/tri.hic.hap1.p_ctg.fa $H/ERR18010482_1.fastq.gz $H/ERR18010482_2.fastq.gz $O/hap1 $L $O/hifi.meryl 48 > hap1.log 2>&1 &
P1=$!
bash $G/code/wgs/scaffold_qc.sh $G/asm25/tri.hic.hap2.p_ctg.fa $H/ERR18010482_1.fastq.gz $H/ERR18010482_2.fastq.gz $O/hap2 $L $O/hifi.meryl 48 > hap2.log 2>&1 &
P2=$!
wait $P1
wait $P2
step "deposited mSicTri.1: compleasm + merqury"
mkdir -p deposited
cd deposited || exit 1
D=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
compleasm run -a $D -o compleasm -l mammalia --odb odb10 -L $L -t 64 > compleasm.log 2>&1
cat compleasm/summary.txt
mkdir -p merqury
cd merqury || exit 1
merqury.sh $O/hifi.meryl $D mq > merqury.log 2>&1
cat mq.qv mq.completeness.stats
step done
