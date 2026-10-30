#!/bin/bash
set -u
cd /home/toki/gof || exit 1
export PATH=/home/toki/toga2/env_mlc2/bin:$PATH
cat /home/toki/Sicista2026/reference/assembly_metadata.txt; grep -v "^## *$" /home/toki/Sicista2026/reference/assembly_stats.txt | head -40
twoBitToFa /home/toki/toga2/inputs/trizona_nodot.2bit stdout | python3 -c '
import sys, re, collections
seqs = {}; name = None; buf = []
for l in sys.stdin:
    if l.startswith(">"):
        if name: seqs[name] = "".join(buf)
        name = l[1:].strip(); buf = []
    else: buf.append(l.strip())
seqs[name] = "".join(buf)
lens = collections.Counter(); per = collections.Counter(); total = 0; near_end = 0
for n, s in seqs.items():
    for m in re.finditer("[Nn]+", s):
        L = len(m.group()); lens[L] += 1; per[n] += 1; total += L
print("contigs", len(seqs), "N runs", sum(lens.values()), "N bases", total)
print("run length counts (top):", lens.most_common(8))
print("runs per chromosome (top):", per.most_common(5))
'
