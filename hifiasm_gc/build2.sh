#!/bin/bash
# v2 = v1 + HIFIASM_GC_DEBUG instrumentation (inactive unless the env var is set).
set -euo pipefail
D=/home/toki/gof/hifiasm_gc
E=/home/toki/miniforge3/envs/asm25
rm -rf $D/gc2; mkdir -p $D/gc2
cd $D/gc2 || exit 1
tar xzf $D/hifiasm_ec9a8b2.tar.gz
patch -p1 < $D/hifiasm_gc_corr_v2.diff
cmp Overlaps.cpp $D/Overlaps.gc2.cpp
make -j 16 CXX=$E/bin/x86_64-conda-linux-gnu-g++ CC=$E/bin/x86_64-conda-linux-gnu-gcc CPPFLAGS="-include pthread.h" INCLUDES="-I$E/include" LIBS="-L$E/lib -Wl,-rpath,$E/lib -lz -lpthread -lm" > make.log 2>&1 || { grep -E " error" make.log | head; exit 1; }
ls -la $D/gc2/hifiasm
# toy: v2 off must equal stock; v2 on must equal v1 on
T=$D/toy
cd $T || exit 1
export TMPDIR=/home/toki/tmp
run() { local n=$1 b=$2 e=$3
  rm -rf $T/$n; mkdir -p $T/$n; cd $T/$n || exit 1
  env HIFIASM_GC_CORR=$e $b -o $T/$n/t -t 16 $T/in/hifi.fa > hifi.log 2>&1
  env HIFIASM_GC_CORR=$e $b -o $T/$n/t -t 16 --h1 $T/in/h1.fq.gz --h2 $T/in/h2.fq.gz $T/in/hifi.fa > hic.log 2>&1
  cd $T || exit 1; }
run off2 $D/gc2/hifiasm off
run on2 $D/gc2/hifiasm on
for p in stock:off2 on:on2; do a=${p%%:*}; b=${p#*:}; n=0; d=0
  for f in $(cd $a && ls *.gfa *.bed); do n=$((n+1)); cmp -s $a/$f $b/$f || { d=$((d+1)); echo "DIFF $a $b $f"; }; done
  echo "$a vs $b: $n files, $d differ"; done
echo build2_done
