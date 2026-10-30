#!/bin/bash
# After wgs_hifiasm.sh: fill trizona N runs from hap1/hap2 contigs, map all reads to each
# haplotype, count junction reads. Waits for the hifiasm wrapper PID given as $1.
set -euo pipefail
W=/home/toki/gof/wgs
cd $W || exit 1
[ "$(pwd)" = "$W" ] || exit 1
export PATH=/home/toki/toga2/env_mlc2/bin:/home/toki/miniforge3/envs/sinederella/bin:$PATH
export TMPDIR=/home/toki/tmp
C=/home/toki/gof/code/wgs
QFA=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
READS=/home/toki/gof/reads/ERR16751478.fastq.gz
step() { echo "#### $* -- $(date '+%F %T')"; }
step "waiting for hifiasm wrapper $1"
while kill -0 "$1" 2> /dev/null; do sleep 300; done
grep -q '^done' run.log || { echo "hifiasm wrapper did not finish (run.log)"; exit 1; }
for h in hap1 hap2; do [ -s tri.bp.$h.p_ctg.fa ] || { echo "missing tri.bp.$h.p_ctg.fa"; exit 1; }; done
step gapfill
python3 $C/gapfill.py $QFA gapfill hap1=tri.bp.hap1.p_ctg.fa hap2=tri.bp.hap2.p_ctg.fa --threads 32
cut -f6 gapfill/gaps.tsv | sort | uniq -c
for h in hap1 hap2; do
  step "map reads to $h"
  minimap2 -x map-hifi -t 64 tri.bp.$h.p_ctg.fa $READS 2> map_$h.log > reads_vs_$h.paf
done
step readcheck
python3 $C/readcheck.py gapfill hap1=reads_vs_hap1.paf hap2=reads_vs_hap2.paf
step finished
