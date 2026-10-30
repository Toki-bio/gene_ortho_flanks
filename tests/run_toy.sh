#!/bin/bash
set -euo pipefail
C=/home/toki/gof/code
W=/home/toki/gof/toy
mkdir -p $W
cd $W || exit 1
export PATH=/home/toki/toga2/env_mlc2/bin:/home/toki/miniforge3/envs/sinederella/bin:$PATH
export TMPDIR=/home/toki/tmp
T=/home/toki/toga2/b2t
TO=$T/full/toga_out
python3 $C/tests/make_toy.py $TO $T/sbe_ref_isoforms.tsv /home/toki/toga2/inputs/trizona_nodot.2bit $W 100
for v in keep del gap; do
  faToTwoBit $v.fa $v.2bit
  minimap2 -x asm20 -d $v.mmi $v.fa 2> $v.idx.log
  python3 $C/gof.py --toga $TO --ref-bed $T/sbe_ref_annotation.bed --isoforms $T/sbe_ref_isoforms.tsv \
    --ref-2bit $T/sbe.2bit --query-2bit $v.2bit --query-index $v.mmi --genes genes.txt --out $W/out_$v \
    --threads 8 --no-neighbours
  echo "== $v: f1000 / f10000 verdicts"
  awk -F'\t' 'NR>1{print $5"\t"$8}' out_$v/verdicts.tsv | sort | uniq -c
done
