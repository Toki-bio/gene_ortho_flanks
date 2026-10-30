#!/bin/bash
# After run_toy.sh (PID given) ends: rerun gof.py on the toy versions with the existing WGA PAFs, then the real run.
set -u
while kill -0 "$1" 2>/dev/null; do sleep 30; done
export PATH=/home/toki/toga2/env_mlc2/bin:/home/toki/miniforge3/envs/sinederella/bin:$PATH
C=/home/toki/gof/code; T=/home/toki/toga2/b2t; TO=$T/full/toga_out
cd /home/toki/gof/toy || exit 1
for v in keep del gap; do
  python3 $C/gof.py --toga $TO --ref-bed $T/sbe_ref_annotation.bed --isoforms $T/sbe_ref_isoforms.tsv \
    --ref-2bit $T/sbe.2bit --query-2bit $v.2bit --query-index $v.mmi --genes genes.txt --out out_$v \
    --threads 8 --no-neighbours --wga-paf wga_$v.paf > /dev/null
  echo "== toy $v: gene call"; awk -F'\t' 'NR>1{print $NF}' out_$v/verdicts.tsv | sort | uniq -c
done
bash $C/tests/run_b2t.sh > /home/toki/gof/b2t.log 2>&1
sed -n '/== gene call/,$p' /home/toki/gof/b2t.log
bash $C/tests/empties.sh
