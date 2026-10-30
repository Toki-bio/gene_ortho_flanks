#!/bin/bash
# 20 kb on each side of the trizona 100-N gaps that face the betulina "missing" stretches:
# GC, self-similar (tandem/duplicated) bases, longest exact short-period repeat.
cd /home/toki/gof/b2t || exit 1
export PATH=/home/toki/toga2/env_mlc2/bin:/home/toki/miniforge3/envs/sinederella/bin:$PATH
Q=/home/toki/toga2/inputs/trizona_nodot.2bit
awk -F'\t' 'NR>1 && $10=="deletion_larger_than_flanks" && $12>0 {print $11}' no_anchor_profile_n.tsv | sort -u > gap_intervals.txt
while read loc; do
  c=${loc%%:*}; se=${loc#*:}; s=${se%-*}; e=${se#*-}
  twoBitToFa $Q:$c:$s-$e stdout | python3 -c '
import sys, re
seq = "".join(sys.stdin.read().split("\n")[1:])
for m in re.finditer("N+", seq): print(m.start(), m.end())' | while read a b; do
    gs=$((s + a)); ge=$((s + b))
    for side in L R; do
      if [ $side = L ]; then x=$((gs - 20000)); y=$gs; else x=$ge; y=$((ge + 20000)); fi
      twoBitToFa $Q:$c:$x-$y /home/toki/tmp/edge.fa
      dup=$(minimap2 -c -x asm20 -X -N 50 -p 0.1 /home/toki/tmp/edge.fa /home/toki/tmp/edge.fa 2>/dev/null | awk '!($3==$8&&$4==$9) && $11>=500 && $10/$11>=0.9 {b+=$11} END{print b+0}')
      python3 - "$c:$gs" "$side" "$dup" <<'PY'
import re, sys
seq = "".join(open("/home/toki/tmp/edge.fa").read().split("\n")[1:]).upper()
acgt = sum(seq.count(b) for b in "ACGT") or 1
gc = (seq.count("G") + seq.count("C")) / acgt
runs = [len(m.group()) for m in re.finditer(r"((\w{1,6})\2{4,})", seq)]
print(f"{sys.argv[1]}\t{sys.argv[2]}\tGC {gc:.3f}\tself-similar bp {sys.argv[3]}\tshort-period repeat bp {sum(runs)}\tlongest {max(runs, default=0)}")
PY
    done
  done
done < gap_intervals.txt
rm -f /home/toki/tmp/edge.fa
