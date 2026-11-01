#!/bin/bash
set -euo pipefail
D=/home/toki/gof/hifiasm_gc
cd $D || exit 1
for v in stock gc; do
  rm -rf $D/$v; mkdir -p $D/$v
  cd $D/$v || exit 1
  tar xzf $D/hifiasm_ec9a8b2.tar.gz
  if [ $v = gc ]; then patch -p1 < $D/hifiasm_gc_corr.diff; cmp Overlaps.cpp $D/Overlaps.gc.cpp; fi
  E=/home/toki/miniforge3/envs/asm25
  make -j 16 CXX=$E/bin/x86_64-conda-linux-gnu-g++ CC=$E/bin/x86_64-conda-linux-gnu-gcc CPPFLAGS="-include pthread.h" INCLUDES="-I$E/include" LIBS="-L$E/lib -Wl,-rpath,$E/lib -lz -lpthread -lm" > make.log 2>&1 || { tail -30 make.log; exit 1; }
  ./hifiasm --version
  cd $D || exit 1
done
ls -la $D/stock/hifiasm $D/gc/hifiasm
/home/toki/miniforge3/envs/asm25/bin/hifiasm --version
