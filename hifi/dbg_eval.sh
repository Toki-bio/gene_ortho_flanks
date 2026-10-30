#!/bin/bash
cd /home/toki/gof/hifi_toy || exit 1
export PATH=/home/toki/miniforge3/envs/sinederella/bin:/home/toki/toga2/env_mlc2/bin:$PATH
python3 - <<'PY'
seqs={};n=None
for l in open("bait.fa"):
    l=l.strip()
    if l.startswith(">"): n=l[1:]; seqs[n]=[]
    else: seqs[n].append(l)
s="".join(seqs["R00|q"]); print("R00|q len", len(s), "N at", s.upper().find("N"), "count", s.upper().count("N"))
open("/home/toki/tmp/r00q.fa","w").write(">R00q\n"+s+"\n")
PY
grep ">" asm/R00/contigs.fa | head; for f in asm/R00/contigs.fa; do awk '/^>/{n=$1} !/^>/{print n, length($0)}' $f; done
minimap2 -c -x asm5 --secondary=no /home/toki/tmp/r00q.fa asm/R00/contigs.fa 2>/dev/null | cut -f1-11
