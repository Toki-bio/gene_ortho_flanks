#!/bin/bash
# Whole-genome hifiasm of S. trizona HiFi (ERR16751478), no purging (-l0), both haplotypes kept.
set -euo pipefail
W=/home/toki/gof/wgs
cd $W || exit 1
[ "$(pwd)" = "$W" ] || exit 1
export PATH=/home/toki/miniforge3/envs/hifiasm/bin:$PATH
export TMPDIR=/home/toki/tmp
echo "start $(date '+%F %T') $(hifiasm --version)"
hifiasm -o $W/tri -t 64 -l0 /home/toki/gof/reads/ERR16751478.fastq.gz > $W/hifiasm.log 2>&1
for p in bp.p_ctg bp.hap1.p_ctg bp.hap2.p_ctg; do
  awk '$1=="S"{print ">"$2"\n"$3}' $W/tri.$p.gfa > $W/tri.$p.fa
done
echo "done $(date '+%F %T')"
