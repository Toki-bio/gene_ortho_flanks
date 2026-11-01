#!/bin/bash
# Genome-wide independent assembly of the S. trizona HiFi reads with HiCanu (canu 2.3 -pacbio-hifi, bioconda env
# /home/toki/miniforge3/envs/canu). Local execution, <= 56 threads (shared server). A disk watchdog stops the run
# (own PID tree only) if /home free space drops below 1.1 TB.
# Usage: asm2_hicanu.sh   (run with setsid nohup; log /home/toki/gof/asm2/hicanu.log)
set -euo pipefail
cd /home/toki/gof/asm2 || exit 1
export TMPDIR=/home/toki/tmp PATH=/home/toki/miniforge3/envs/canu/bin:$PATH
R=/home/toki/gof/reads/ERR16751478.fastq.gz
echo "#### start $(date '+%F %T') pid $$"
df -h /home | tail -n 1
(
  while sleep 300; do
    free_kb=$(df -k --output=avail /home | tail -n 1)
    if [ "$free_kb" -lt 1181116006 ]; then              # 1.1 TiB
      echo "#### WATCHDOG $(date '+%F %T'): /home free ${free_kb} kB < 1.1 TB, stopping canu"
      pkill -TERM -s $$ || true                         # everything in this script's own session (setsid)
      exit 0
    fi
  done
) &
WD=$!
/usr/bin/time -v canu -p tri -d hicanu genomeSize=2.8g useGrid=false maxThreads=56 maxMemory=240g \
  -pacbio-hifi $R > hicanu.out 2> hicanu.time || echo "canu exit $?"
kill $WD 2>/dev/null || true
grep -E "Elapsed|Maximum res|Exit" hicanu.time || true
ls -la hicanu/ | head -30
if [ -s hicanu/tri.contigs.fasta ]; then
  awk '/^>/{n++} !/^>/{s+=length($0)} END{print n" contigs "s" bp"}' hicanu/tri.contigs.fasta
  /home/toki/miniforge3/envs/asm25/bin/gfastats hicanu/tri.contigs.fasta 2>/dev/null | head -n 20 || true
fi
df -h /home | tail -n 1
echo "#### done $(date '+%F %T')"
