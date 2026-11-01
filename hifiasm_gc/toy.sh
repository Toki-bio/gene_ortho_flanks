#!/bin/bash
# Toy: fished reads of the 10 studied regions + 4 controls (ends/sv/fished.fa) and the first 1M Hi-C pairs.
# Runs stock, patched (off) and patched (on); compares the output files.
set -euo pipefail
D=/home/toki/gof/hifiasm_gc
T=$D/toy
mkdir -p $T/in
cd $T || exit 1
export TMPDIR=/home/toki/tmp
if [ ! -s in/h2.fq.gz ]; then
  cp /home/toki/gof/ends/sv/fished.fa in/hifi.fa
  for r in 1 2; do { zcat /home/toki/gof/hic/ERR18010482_$r.fastq.gz || true; } | head -4000000 | gzip -1 > in/h$r.fq.gz; done
fi
run() { # name binary env
  local n=$1 b=$2 e=$3
  rm -rf $T/$n; mkdir -p $T/$n; cd $T/$n || exit 1
  env HIFIASM_GC_CORR=$e $b -o $T/$n/t -t 16 $T/in/hifi.fa > hifi.log 2>&1
  env HIFIASM_GC_CORR=$e $b -o $T/$n/t -t 16 --h1 $T/in/h1.fq.gz --h2 $T/in/h2.fq.gz $T/in/hifi.fa > hic.log 2>&1
  cd $T || exit 1
}
run stock $D/stock/hifiasm off
run off $D/gc/hifiasm off
run on $D/gc/hifiasm on
for n in off on; do
  echo "== stock vs $n"
  for f in $(cd stock && ls *.gfa *.bed | sort); do
    if cmp -s stock/$f $n/$f; then echo "same $f"; else echo "DIFF $f"; fi
  done
done
grep gc_corr stock/hic.log off/hic.log on/hic.log || true
grep "primary contig coverage" on/hic.log || true
for n in stock off on; do for h in hap1 hap2; do
  awk -v n=$n -v h=$h '$1=="S"{c++; s+=length($3)} END{print n, h, c, s}' $n/t.hic.$h.p_ctg.gfa
done; done
