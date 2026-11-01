#!/bin/bash
# Map the fished reads to the h1 and alt windows separately, TRF on the h1 windows, then recheck_reads.py.
# Usage: recheck_reads.sh DIR (holds win.fa + fished.fq from recheck_fish.sh) [THREADS]
set -euo pipefail
D=$1
T=${2:-24}
cd "$D" || exit 1
C=/home/toki/gof/code/hifi
export PATH=/home/toki/miniforge3/envs/sinederella/bin:$PATH
awk '/^>/{p=($0 ~ /\|h1\|/)} p' win.fa > h1.fa
awk '/^>/{p=($0 ~ /\|alt\|/)} p' win.fa > alt.fa
echo "#### map $(date '+%F %T')"
minimap2 -c --cs -x map-hifi --secondary=no -t $T h1.fa fished.fq > h1.paf 2> mm_h1.log
if [ -s alt.fa ]; then minimap2 -c --cs -x map-hifi --secondary=no -t $T alt.fa fished.fq > alt.paf 2> mm_alt.log; else : > alt.paf; fi
echo "#### trf $(date '+%F %T')"
trf h1.fa 2 7 7 80 10 50 2000 -h -ngs > h1.trf.dat || true
echo "#### reads $(date '+%F %T')"
python3 $C/recheck_reads.py win.fa h1.paf alt.paf fished.fq h1.trf.dat res
echo "#### done $(date '+%F %T')"
