#!/bin/bash
# Full trizona Hi-C phasing with the GC-corrected hifiasm, reusing asm25's ec/ovlp/Hi-C bins (copied, not linked).
set -euo pipefail
D=/home/toki/gof/hifiasm_gc/asm_gc
A=/home/toki/gof/asm25
mkdir -p $D
cd $D || exit 1
[ "$(pwd)" = "$D" ] || exit 1
export TMPDIR=/home/toki/tmp
echo "#### copy bins -- $(date '+%F %T')"
for f in tri.ec.bin tri.ovlp.source.bin tri.ovlp.reverse.bin tri.hic.lk.bin tri.hic.tlb.bin; do
  if [ ! -s $f ] || [ $(stat -c %s $f) -ne $(stat -c %s $A/$f) ]; then cp $A/$f $D/$f; fi
done
ls -la $D
df -h /home | tail -1
echo "#### hifiasm GC-corrected Hi-C phasing -- $(date '+%F %T')"
HIFIASM_GC_CORR=on /home/toki/gof/hifiasm_gc/gc/hifiasm -o $D/tri -t 48 \
  --h1 /home/toki/gof/hic/ERR18010482_1.fastq.gz --h2 /home/toki/gof/hic/ERR18010482_2.fastq.gz \
  /home/toki/gof/reads/ERR16751478.fastq.gz > hic.log 2>&1
for p in hic.hap1.p_ctg hic.hap2.p_ctg; do
  awk '$1=="S"{print ">"$2"\n"$3}' $D/tri.$p.gfa > $D/tri.$p.fa
done
grep -c '>' tri.hic.hap1.p_ctg.fa tri.hic.hap2.p_ctg.fa
echo "#### done -- $(date '+%F %T')"
