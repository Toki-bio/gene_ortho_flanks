#!/usr/bin/env python3
"""Depth obs/exp by kb from the contig end, median per cause, for ends of contigs >= MINLEN; plus where hemizygous
contigs sit in trizona (share of end flanks mapping to the trizona X).
Usage: end_ramp.py END_CAUSES.tsv ENDS.bedcov CURVE [MINLEN=1000000] [TRI_X=OZ418353.1]"""
import statistics
import sys
from collections import Counter, defaultdict

from calib import load_curve

ec, bcov, curve = sys.argv[1:4]
MINLEN = int(sys.argv[4]) if len(sys.argv) > 4 else 1000000
TX = sys.argv[5] if len(sys.argv) > 5 else "OZ418353.1"
expected, REF = load_curve(curve)
rows = [l.rstrip("\n").split("\t") for l in open(ec)]
h = {k: i for i, k in enumerate(rows[0])}
info = {r[0]: r for r in rows[1:]}
dep = defaultdict(dict)
for l in open(bcov):
    c, s, e, n, v = l.split()
    ctg, side, k = n.split("|")
    dep[f"{ctg}|{side}"][int(k)] = float(v) / (int(e) - int(s))
KS = [0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 25, 29]
by = defaultdict(lambda: defaultdict(list))
for e, r in info.items():
    if int(r[h['ctg_len']]) < MINLEN:
        continue
    ex = expected(round(float(r[h['gc10k']]))) * REF
    for k in KS:
        if k in dep[e]:
            by[r[h['cause']]][k].append(dep[e][k] / ex)
print("cause\tn\t" + "\t".join(f"kb{k}" for k in KS))
for c, d in sorted(by.items(), key=lambda x: -len(x[1][0])):
    print(f"{c}\t{len(d[0])}\t" + "\t".join(f"{statistics.median(d[k]):.2f}" for k in KS))
hemi = Counter()
for e, r in info.items():
    if r[h['hemi']] == "1":
        L = int(r[h['ctg_len']])
        size = ">=1Mb" if L >= 1e6 else "<1Mb"
        hemi[(size, "triX" if r[h['tri_chrom']] == TX else "unmapped" if not r[h['tri_chrom']] else "autosome")] += 1
print("hemizygous-contig ends by trizona chromosome:", dict(sorted(hemi.items())))
