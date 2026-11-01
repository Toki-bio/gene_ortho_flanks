#!/bin/bash
# Negative control for overlap merges: the two toy gaps accepted with the read-measured E are re-planned with E off by
# -2 kb and +-500 bp; with the reads' span rule these should be rejected (or show which tolerance applies).
set -euo pipefail
cd /home/toki/gof/asm2/ovtoy || exit 1
C=/home/toki/gof/code/hifi
D1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
export PATH=/home/toki/miniforge3/envs/hifiasm/bin:/home/toki/miniforge3/envs/jf/bin:/home/toki/miniforge3/envs/lrtools/bin:$PATH
for d in -2000 -500 500; do
  O=neg$d
  python3 $C/asm2_overlap.py $D1 /home/toki/gof/gf_all/dep1/gaps.tsv $O "gap0112:$((14945+d)):536:17439,gap0179:$((5276+d)):4689:1239" > /dev/null 2>&1
  n=$(grep -c '>' $O/cand.fa)
  jellyfish count -m 25 -C -s $((4 * n + 1000)) -t 8 --if $O/cand.fa -o $O/cand.jf $D1
  jellyfish dump -c $O/cand.jf > $O/cand.counts
  python3 $C/repair_gaps.py words $O $O/cand.counts 2>/dev/null
  python3 $C/repair_gaps.py verify $D1 $O fished.fq dep1 8 --edited-only > /dev/null
  rm -f $O/*.fa.tmp $O/dep1.repaired.fa $O/dep1.genesearch.fa $O/cand.jf $O/cand.counts
  echo "== E offset $d"; cut -f2,10,12-14 $O/repairs.tsv | tail -n +2
done
