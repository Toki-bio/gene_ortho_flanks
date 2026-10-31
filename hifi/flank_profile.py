#!/usr/bin/env python3
"""Read depth along the flanks of each region, split by alignment identity, to see where the reads of the copy
WITHOUT the piece stop. The flanks come from the copy WITH the piece, so its reads align at ~HiFi accuracy and reads
of the other copy at about 0.3 % lower (hap1 vs hap2 divergence measured 2026-10-07: 3.03 SNVs/kb).
Input: OUT_DIR of classify_ends.py (parts.fa) + the same reads. Prints per region and flank, in 1 kb windows counted
from the gap: reads with identity >= HI, reads below HI (each read once per window, primary chains only).
Usage: flank_profile.py CLASSIFY_OUT_DIR READS.fq [THREADS] [HI]"""
import os
import subprocess
import sys
from collections import defaultdict

d, reads = sys.argv[1:3]
threads = sys.argv[3] if len(sys.argv) > 3 else "8"
HI = float(sys.argv[4]) if len(sys.argv) > 4 else 0.998
W = 1000
paf = os.path.join(d, "parts_vs_reads.paf")
if not os.path.exists(paf):
    subprocess.run(f"minimap2 -c -x map-hifi --secondary=no -t {threads} {d}/parts.fa {reads} > {paf} 2> /dev/null",
                   shell=True, check=True)
tlen = {}
for l in open(paf):
    f = l.split("\t")
    tlen[f[5]] = int(f[6])
cov = defaultdict(lambda: defaultdict(lambda: [set(), set()]))
for l in open(paf):
    f = l.split("\t")
    t, ts, te, m, al = f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])
    if not (t.endswith("|L") or t.endswith("|R")) or al < 2000:
        continue
    hi = m / al >= HI
    for w in range(ts // W, (te - 1) // W + 1):
        if w * W >= ts and (w + 1) * W <= te:          # read covers the whole window
            dist = tlen[t] - (w + 1) * W if t.endswith("|L") else w * W   # window distance from the gap
            cov[t][dist // W][0 if hi else 1].add(f[0])
for t in sorted(cov):
    ks = sorted(cov[t])
    print(f"{t}  kb_from_gap: " + " ".join(f"{k}:{len(cov[t][k][0])}/{len(cov[t][k][1])}" for k in ks if k < 10))
