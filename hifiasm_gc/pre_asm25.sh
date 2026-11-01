#!/bin/bash
# asm25-side parts of eval.sh run early (same output names, eval.sh skips them when present).
set -euo pipefail
G=/home/toki/gof
O=$G/hifiasm_gc/eval
A=$G/asm25
mkdir -p $O
cd $O || exit 1
export PATH=/home/toki/miniforge3/envs/asm25/bin:$PATH
export MERQURY=$(ls -d /home/toki/miniforge3/envs/asm25/share/merqury*)
export TMPDIR=/home/toki/tmp
T=16
for h in hap1 hap2; do
  [ -s cpl_asm25_$h/summary.txt ] || compleasm run -a $A/tri.hic.$h.p_ctg.fa -o cpl_asm25_$h -l mammalia --odb odb10 -L $G/compleasm_lib -t $T > cpl_asm25_$h.log 2>&1
  grep -E '^[SDFM]:' cpl_asm25_$h/summary.txt
done
mkdir -p mq_asm25; cd $O/mq_asm25 || exit 1
ln -sf $A/tri.hic.hap1.p_ctg.fa asm25_h1.fa; ln -sf $A/tri.hic.hap2.p_ctg.fa asm25_h2.fa
[ -s mq.qv ] || OMP_NUM_THREADS=$T merqury.sh $G/sqc/hifi.meryl asm25_h1.fa asm25_h2.fa mq > merqury.log 2>&1
cat mq.completeness.stats mq.qv
echo pre_done
