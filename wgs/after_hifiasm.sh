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
# hifiasm 0.16.1 without Hi-C/trio writes primary (p_ctg) + alternate (a_ctg), not hap1/hap2;
# with -l0 the primary set keeps both alleles of most heterozygous regions.
for h in p a; do
  [ -s tri.bp.${h}_ctg.fa ] || awk '$1=="S"{print ">"$2"\n"$3}' tri.bp.${h}_ctg.gfa > tri.bp.${h}_ctg.fa
  [ -s tri.bp.${h}_ctg.fa ] || { echo "missing tri.bp.${h}_ctg.fa"; exit 1; }
done
step gapfill
python3 $C/gapfill.py $QFA gapfill pri=tri.bp.p_ctg.fa alt=tri.bp.a_ctg.fa --threads 32
cut -f6 gapfill/gaps.tsv | sort | uniq -c
for h in p a; do
  step "map reads to ${h}_ctg"
  minimap2 -x map-hifi -t 64 tri.bp.${h}_ctg.fa $READS 2> map_$h.log > reads_vs_$h.paf
done
step readcheck
python3 $C/readcheck.py gapfill pri=reads_vs_p.paf alt=reads_vs_a.paf
step finished
