#!/bin/bash
# Self-alignment of each betulina stretch facing a trizona assembly gap: bases in
# internal duplicates (>= 2 kb, >= 95 % identity), longest duplicate.
cd /home/toki/gof/b2t
export PATH=/home/toki/toga2/env_mlc2/bin:/home/toki/miniforge3/envs/sinederella/bin:$PATH
for x in ptg000158l:6649030-6675945 ptg000006l:3730556-4015233 ptg000045l:499662-538368 ptg000153l:15200357-15285276 \
         ptg000116l:33114307-33193032 ptg000034l:6757457-6784600 ptg000116l:31889650-31996206 ptg000365l:10912552-11001798 ptg000034l:18860251-18886176; do
  twoBitToFa /home/toki/toga2/b2t/sbe.2bit:$x /home/toki/tmp/self.fa
  minimap2 -c -x asm20 -X -N 50 -p 0.1 /home/toki/tmp/self.fa /home/toki/tmp/self.fa 2>/dev/null \
    | awk -v x=$x 'BEGIN{OFS="\t"} $1==$6 && !($3==$8 && $4==$9) && $11>=2000 && $10/$11>=0.95 {n++; b+=$11; if($11>m)m=$11} END{print x, "dup_pairs", n+0, "dup_bp", b+0, "longest", m+0}'
done
rm -f /home/toki/tmp/self.fa
