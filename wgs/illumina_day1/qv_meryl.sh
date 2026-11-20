#!/bin/bash
# Assembly QV from Illumina k-mers, Merqury formula, with meryl only (merqury.sh is not installed).
# usage: qv_meryl.sh asm.fa R1.fq.gz R2.fq.gz tag [threads]
# reads.meryl is built once and reused. QV = -10 log10(1 - (1 - missing/total)^(1/k)), where total = all k-mer positions in the
# assembly and missing = assembly k-mer positions whose k-mer never occurs in the reads.
set -euo pipefail
A=$1; R1=$2; R2=$3; TAG=$4; T=${5:-32}; K=21
MER=/usr/local/bin/meryl
[ -d reads.meryl ] || $MER count k=$K threads=$T memory=200 $R1 $R2 output reads.meryl > /dev/null 2>&1
rm -rf asm_$TAG.meryl miss_$TAG.meryl
$MER count k=$K threads=$T memory=100 $A output asm_$TAG.meryl > /dev/null 2>&1
$MER difference asm_$TAG.meryl reads.meryl output miss_$TAG.meryl > /dev/null 2>&1
tot=$($MER statistics asm_$TAG.meryl 2>/dev/null | awk '$1=="present"{print $2; exit}')
mis=$($MER statistics miss_$TAG.meryl 2>/dev/null | awk '$1=="present"{print $2; exit}')
awk -v t=$tot -v m=$mis -v k=$K -v tag=$TAG 'BEGIN{e=1-(1-m/t)^(1/k); printf "%s: assembly k-mers %d, missing from reads %d, QV %.2f (error rate %.2e)\n", tag, t, m, (e>0?-10*log(e)/log(10):99), e}'
rm -rf asm_$TAG.meryl miss_$TAG.meryl
