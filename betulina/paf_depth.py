#!/usr/bin/env python3
"""Mean read depth per bin from a PAF on stdin (primary alignments, MAPQ >= Q), like mosdepth --by BIN.
Output (BED, every bin of every sequence in FAI): contig  start  end  mean_depth
Usage: minimap2 ... | paf_depth.py GENOME.fa.fai BIN [Q] > depth.bed"""
import sys
from collections import defaultdict

fai, BIN = sys.argv[1], int(sys.argv[2])
Q = int(sys.argv[3]) if len(sys.argv) > 3 else 20
ln = {l.split()[0]: int(l.split()[1]) for l in open(fai)}
cov = defaultdict(int)
for l in sys.stdin:
    f = l.split("\t")
    if int(f[11]) < Q or "tp:A:P" not in l:
        continue
    t, a, e = f[5], int(f[7]), int(f[8])
    while a < e:
        b = a // BIN * BIN; k = min(e, b + BIN) - a
        cov[(t, b)] += k; a += k
for t, L in ln.items():
    for b in range(0, L, BIN):
        w = min(BIN, L - b)
        print(f"{t}\t{b}\t{b + w}\t{cov[(t, b)] / w:.2f}")
