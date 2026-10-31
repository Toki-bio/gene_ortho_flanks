#!/usr/bin/env python3
"""k-mer copy number of the gap flanks per kb of distance from the gap, median per cause group (GC-calibrated per flank).
A dropout on one copy at the gap predicts a ramp: ~1.0 far away, ~0.5 next to the gap (molecules spanning the
unreadable spot give no reads, so that copy's coverage fades over about one read length).
Usage: flank_ramp.py FLANKS.fa FLANKS.wig CALIB GAP_CAUSES.tsv [PEAK]"""
import statistics
import sys
from collections import defaultdict

fa, wig, calib, causes = sys.argv[1:5]
PEAK = float(sys.argv[5]) if len(sys.argv) > 5 else 26
exp = {int(g): float(v) for g, v in (x.split(":") for x in open(calib).readline().split(":", 1)[1].split())}
expected = lambda gc: exp[min(exp, key=lambda k: abs(k - gc))]
seqs, name = {}, None
for l in open(fa):
    if l.startswith(">"):
        name = l[1:].strip(); seqs[name] = ""
    else:
        seqs[name] += l.strip().upper()
cause = {l.split("\t")[0]: l.rstrip("\n").split("\t")[-1] for l in list(open(causes))[1:]}
vals = defaultdict(lambda: defaultdict(list))
cur = None
for l in open(wig):
    if l.startswith("variableStep"):
        cur = l.split("chrom=")[1].split()[0]; continue
    p, c = map(int, l.split())
    if c > 0:
        vals[cur][p // 1000 if cur.endswith("|R") else None].append(c) if cur.endswith("|R") else \
            vals[cur][(len(seqs[cur]) - 1 - p) // 1000].append(c)
ramp = defaultdict(lambda: defaultdict(list))
for n, d in vals.items():
    s = seqs[n]
    gc = 100 * (s.count("G") + s.count("C")) / max(1, sum(s.count(x) for x in "ACGT"))
    e = expected(round(gc))
    g = cause.get(n.split("|")[0], "?")
    for k, v in d.items():
        if k is not None and len(v) >= 250:
            ramp[g][k].append(statistics.median(v) / PEAK / e)
for g in sorted(ramp, key=lambda x: -len(ramp[x][0])):
    print(f"{g:22s} n={len(ramp[g][0]):4d}  obs/exp by kb from gap: " +
          " ".join(f"{k}:{statistics.median(ramp[g][k]):.2f}" for k in sorted(ramp[g]) if k < 10))
