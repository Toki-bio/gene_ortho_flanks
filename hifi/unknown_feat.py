#!/usr/bin/env python3
"""Summarise per sequence group (UNK gap flanks, ONE / TWO copy fills, CTRL segments): repeat classes (RepeatMasker
.out with the rodent Dfam subset), tandem repeats (TRF .dat), GC, k-mer copy number (meryl-lookup wig; median count /
PEAK) and the fraction of k-mers at >= 4x the two-copy peak (collapsed high-copy repeat).
Also per-sequence table OUT/feat_per_seq.tsv.
Usage: unknown_feat.py FEAT.fa RM.out TRF.dat WIG OUT_DIR [PEAK]"""
import os
import statistics
import sys
from collections import defaultdict

fa, rmo, trf, wig, out = sys.argv[1:6]
PEAK = float(sys.argv[6]) if len(sys.argv) > 6 else 26
seqs, name = {}, None
for l in open(fa):
    if l.startswith(">"):
        name = l[1:].split()[0]; seqs[name] = []
    else:
        seqs[name].append(l.strip().upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
grp = lambda n: n.split("_")[0]
cls = defaultdict(lambda: defaultdict(int))
for l in open(rmo):
    f = l.split()
    if len(f) < 15 or not f[0].isdigit():
        continue
    q, a, b, c = f[4], int(f[5]), int(f[6]), f[10]
    cls[q][c.split("/")[0] if c.split("/")[0] in ("LINE", "SINE", "LTR", "DNA") else c] += b - a + 1
tr = defaultdict(int); per = defaultdict(list)
cur = None
for l in open(trf):
    if l.startswith("@"):                      # trf -ngs: "@name" then one line per repeat
        cur = l[1:].split()[0]
    elif cur and l[:1].isdigit():
        f = l.split()
        tr[cur] += int(f[1]) - int(f[0]) + 1
        per[cur].append((int(f[1]) - int(f[0]) + 1, int(f[2]), f[13] if len(f) > 13 else ""))
cnt = defaultdict(list)
for l in open(wig):
    if l.startswith("variableStep"):
        cur = l.split("chrom=")[1].split()[0]; continue
    cnt[cur].append(int(l.split()[1]))
rows = []
with open(os.path.join(out, "feat_per_seq.tsv"), "w") as fh:
    fh.write("seq\tlen\tGC\tkmer_cn\thighcopy_frac\ttandem_frac\ttop_period\trepeat_classes\n")
    for n, s in seqs.items():
        L = len(s)
        gc = 100 * (s.count("G") + s.count("C")) / max(1, sum(s.count(x) for x in "ACGT"))
        v = cnt.get(n, [])
        nz = [x for x in v if x > 0]
        kc = statistics.median(nz) / PEAK if nz else 0
        hi = sum(1 for x in v if x >= 4 * PEAK) / max(1, len(v))
        tf = min(1.0, tr[n] / L)
        tp = max(per[n])[1] if per[n] else 0
        rc = sorted(cls[n].items(), key=lambda kv: -kv[1])
        rcs = ",".join(f"{k}:{100 * x / L:.0f}%" for k, x in rc[:4])
        rows.append((grp(n), n, L, gc, kc, hi, tf, tp, cls[n]))
        fh.write(f"{n}\t{L}\t{gc:.1f}\t{kc:.2f}\t{hi:.3f}\t{tf:.3f}\t{tp}\t{rcs}\n")
print("group  n   GC   kmerCN  highcopy>4x  tandem%  masked% by class (mean over sequences)")
for g in ("UNK", "ONE", "TWO", "CTRL"):
    rr = [r for r in rows if r[0] == g]
    if not rr:
        continue
    tot = defaultdict(float)
    for r in rr:
        for k, x in r[8].items():
            tot[k] += 100 * x / r[2] / len(rr)
    top = " ".join(f"{k}:{v:.1f}" for k, v in sorted(tot.items(), key=lambda kv: -kv[1])[:7])
    print(f"{g:5s} {len(rr):4d} {statistics.median(r[3] for r in rr):5.1f} {statistics.median(r[4] for r in rr):6.2f} "
          f"{statistics.mean(r[5] for r in rr):10.3f} {100 * statistics.mean(r[6] for r in rr):8.1f}  {top}")
