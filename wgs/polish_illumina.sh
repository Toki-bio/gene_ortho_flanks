#!/bin/bash
# Short-read polishing of a long-read assembly with pypolca --careful (changes a base only when the reads
# show no support for the assembly's base, so heterozygous sites are left alone).
# Usage: polish_illumina.sh ASM.fa R1.fq.gz R2.fq.gz OUT_DIR [THREADS]
# Output: OUT_DIR/polished.fa, OUT_DIR/fastp.json, OUT_DIR/pypolca/ (report with substitution/indel counts)
set -euo pipefail
A=$(readlink -f $1); R1=$(readlink -f $2); R2=$(readlink -f $3); O=$4; T=${5:-32}
mkdir -p $O
cd $O || exit 1
export PATH=/home/toki/miniforge3/envs/polish/bin:$PATH
export TMPDIR=/home/toki/tmp
step() { echo "#### $* -- $(date '+%F %T')"; }
step "fastp"
fastp -i $R1 -I $R2 -o trim_1.fq.gz -O trim_2.fq.gz -w 16 -j fastp.json -h fastp.html 2> fastp.log
step "pypolca --careful"
pypolca run -a $A -1 trim_1.fq.gz -2 trim_2.fq.gz -t $T -o pypolca --careful -f > pypolca.log 2>&1
cp pypolca/pypolca_corrected.fasta polished.fa
rm -f trim_1.fq.gz trim_2.fq.gz
cat pypolca/pypolca.report 2>/dev/null || tail -20 pypolca.log
step done
