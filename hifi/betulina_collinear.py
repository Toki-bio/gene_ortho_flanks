#!/usr/bin/env python3
"""Ancestral state along each trizona window from the betulina assembly (S. betulina, sister species ~88-95 % id).
Windows are aligned to the betulina contigs (minimap2 asm20, primary chains >= 2 kb). Prints per window, per 5 kb bin,
the fraction of the bin covered by a betulina alignment and a letter for the betulina contig (same letter = same
contig, '.' = no alignment), then the alignment blocks in window order (window span -> betulina contig:span strand,
identity). A trizona stretch that sits inside one collinear betulina block is ancestral; a stretch that betulina
skips (block ends on both sides, joined on one betulina contig) is an insertion in trizona.
Usage: betulina_collinear.py WINDOWS.fa BETULINA.fa OUT_DIR [THREADS]"""
import os
import string
import subprocess
import sys
from collections import defaultdict

win, bet, out = sys.argv[1:4]
threads = sys.argv[4] if len(sys.argv) > 4 else "16"
BIN = 5000
os.makedirs(out, exist_ok=True)
paf = os.path.join(out, "win_vs_betulina.paf")
if not os.path.exists(paf):
    subprocess.run(f"minimap2 -c -x asm20 --secondary=no -t {threads} {bet} {win} > {paf} 2> /dev/null",
                   shell=True, check=True)
blocks, qlen = defaultdict(list), {}
for l in open(paf):
    f = l.split("\t")
    if int(f[3]) - int(f[2]) < 2000 or "tp:A:P" not in l:
        continue
    qlen[f[0]] = int(f[1])
    blocks[f[0]].append((int(f[2]), int(f[3]), f[5], int(f[7]), int(f[8]), f[4], int(f[9]) / int(f[10])))
for w in sorted(qlen):
    bl = sorted(blocks[w])
    letters, lab = {}, iter(string.ascii_uppercase + string.ascii_lowercase + string.digits)
    for b in bl:
        if b[2] not in letters:
            letters[b[2]] = next(lab, "?")
    cov, who = [], []
    for i in range(qlen[w] // BIN):
        a, z = i * BIN, (i + 1) * BIN
        ov = [(min(z, b[1]) - max(a, b[0]), b[2]) for b in bl if b[0] < z and b[1] > a]
        cov.append(min(1.0, sum(o for o, _ in ov) / BIN))
        who.append(letters[max(ov)[1]] if ov else ".")
    print(f"== {w}  {qlen[w] // 1000} kb")
    print("   contig per 5 kb: " + "".join(who))
    print("   cover  per 5 kb: " + "".join(str(min(9, int(c * 10))) for c in cov))
    for b in bl:
        if b[1] - b[0] >= 10000:
            print(f"   {b[0]:>8}-{b[1]:<8} -> {letters[b[2]]}={b[2]}:{b[3]}-{b[4]} {b[5]}  id {b[6]:.3f}")
