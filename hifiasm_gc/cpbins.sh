#!/bin/bash
set -euo pipefail
D=/home/toki/gof/hifiasm_gc/asm_gc
mkdir -p $D
cd $D || exit 1
for f in tri.ec.bin tri.ovlp.source.bin tri.ovlp.reverse.bin tri.hic.lk.bin tri.hic.tlb.bin; do cp /home/toki/gof/asm25/$f $D/$f; done
ls -la $D; echo copied
