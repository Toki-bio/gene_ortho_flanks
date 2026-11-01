#!/bin/bash
# Toy: HiCanu (canu 2.3 -pacbio-hifi) on the round-1 fished reads (~1,400 junction windows at full local depth).
set -euo pipefail
cd /home/toki/gof/asm2/toy || exit 1
export TMPDIR=/home/toki/tmp PATH=/home/toki/miniforge3/envs/canu/bin:$PATH
/usr/bin/time -v /home/toki/miniforge3/envs/canu/bin/canu -p toy -d hicanu genomeSize=15m useGrid=false \
  maxThreads=32 maxMemory=120g -pacbio-hifi /home/toki/gof/repair/fished.fq.gz > hicanu.out 2> hicanu.time
grep -E "Elapsed|Maximum res|Exit" hicanu.time
ls hicanu
awk '/^>/{n++} !/^>/{s+=length($0)} END{print n" contigs "s" bp"}' hicanu/toy.contigs.fasta
