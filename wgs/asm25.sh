#!/bin/bash
# Phased S. trizona assembly with hifiasm 0.25: HiFi-only first (writes the ec/ovlp bins while the
# Hi-C reads download), then Hi-C mode with the same -o prefix, which reuses the bins and only phases.
# Default purging (phased haplotypes, unlike the -l0 gap-rescue run in wgs_hifiasm.sh).
# Usage: asm25.sh OUT_DIR HIFI.fq.gz HIC_DIR [THREADS]   (HIC_DIR as written by dl_hic.sh)
set -euo pipefail
O=$1; HIFI=$2; H=$3; T=${4:-64}
mkdir -p $O
cd $O || exit 1
[ "$(pwd)" = "$O" ] || exit 1
export PATH=/home/toki/miniforge3/envs/asm25/bin:$PATH
export TMPDIR=/home/toki/tmp
step() { echo "#### $* -- $(date '+%F %T')"; }
step "hifiasm $(hifiasm --version) HiFi only"
hifiasm -o $O/tri -t $T $HIFI > hifi.log 2>&1
step "waiting for Hi-C download in $H"
while [ ! -e $H/download.ok ]; do
  [ -e $H/download.failed ] && { echo "Hi-C download failed"; exit 1; }
  sleep 300
done
step "hifiasm Hi-C phasing"
hifiasm -o $O/tri -t $T --h1 $H/ERR18010482_1.fastq.gz --h2 $H/ERR18010482_2.fastq.gz $HIFI > hic.log 2>&1
for p in hic.hap1.p_ctg hic.hap2.p_ctg hic.p_ctg; do
  awk '$1=="S"{print ">"$2"\n"$3}' $O/tri.$p.gfa > $O/tri.$p.fa
done
grep -c '>' tri.hic.hap1.p_ctg.fa tri.hic.hap2.p_ctg.fa
step done
