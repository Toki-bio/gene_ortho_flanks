#!/bin/bash
# Waiter: when the HiCanu genome run (asm2_hicanu.sh) and the wgs-only round 2 (r2_wgs) have both ended, run repair
# round 2 with HiCanu + hifiasm 0.16.1 -l0 sources -> /home/toki/gof/asm2/r2_hicanu (log r2_hicanu.log).
# Usage: asm2_after_hicanu.sh HICANU_PID R2WGS_PID   (setsid nohup; log /home/toki/gof/asm2/after_hicanu.log)
set -euo pipefail
cd /home/toki/gof/asm2 || exit 1
for p in "$@"; do while kill -0 "$p" 2>/dev/null; do sleep 60; done; done
echo "#### waited $(date '+%F %T')"
A=/home/toki/gof/asm2/hicanu/tri.contigs.fasta
if [ ! -s "$A" ]; then echo "no HiCanu contigs ($A); see hicanu.log"; exit 1; fi
T=32 bash /home/toki/gof/code/hifi/asm2_repair2.sh /home/toki/gof/asm2/r2_hicanu "$A" full > /home/toki/gof/asm2/r2_hicanu.log 2>&1
echo "#### done $(date '+%F %T')"
