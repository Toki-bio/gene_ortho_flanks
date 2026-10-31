#!/usr/bin/env python3
"""GC-matched copy number: is a het-free stretch two identical copies or one copy?

Per 10 kb bin of each window: GC %, k-mer copy number (median meryl count / two-copy peak, from kmer_cn.py's WIG) and
het sites (rule of het_window.py, from its cs PAF). Bins with >= 1 het/kb are certainly two copies; their copy number
per 1 % GC class (median) is the two-copy expectation. Every bin then gets observed / expected: ~1 = two copies,
~0.5 = one copy. Pieces (the gap fills, from FILLS_vs_HAP PAFs and windows.bed) are summarised separately.
Usage: gc_calib.py WINDOWS.fa WINDOWS.bed WIG CS_PAF PAF_hap1 PAF_hap2 [PEAK]"""
import re
import statistics
import sys
from collections import defaultdict

wfa, bed, wig, cspaf, paf1, paf2 = sys.argv[1:7]
PEAK = float(sys.argv[7]) if len(sys.argv) > 7 else 26
BIN, FL = 10000, 10500
seqs, name, buf = {}, None, []
for l in list(open(wfa)) + [">"]:
    if l.startswith(">"):
        if name:
            seqs[name] = "".join(buf).upper()
        name, buf = l[1:].split()[0] if len(l) > 1 else None, []
    else:
        buf.append(l.strip())
cnt = defaultdict(lambda: defaultdict(list))
for l in open(wig):
    if l.startswith("variableStep"):
        chrom = l.split("chrom=")[1].split()[0]; continue
    p, c = l.split()
    if int(c) > 0:
        cnt[chrom][int(p) // BIN].append(int(c))
# het sites per bin
tok = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+)")
cov = defaultdict(lambda: defaultdict(int)); subs = defaultdict(lambda: defaultdict(int))
for l in open(cspaf):
    f = l.rstrip("\n").split("\t")
    t, ts, te = f[5], int(f[7]), int(f[8])
    if te - ts < 2000:
        continue
    for b in range(ts // 1000, te // 1000 + 1):
        cov[t][b] += 1
    cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
    p = ts
    for m in tok.findall(cs):
        if m[0] == ":":
            p += int(m[1:])
        elif m[0] == "*":
            subs[(t, p)][m[2]] += 1; p += 1
        elif m[0] == "-":
            p += len(m) - 1
het = defaultdict(int)
for (t, p), d in subs.items():
    n, c = max(d.values()), cov[t][p // 1000]
    if c >= 6 and n >= 3 and 0.25 <= n / c <= 0.75:
        het[(t, p // BIN)] += 1
# piece spans in window coordinates
wins = {}
for l in open(bed):
    h, c, s, e, lab = l.split()
    wins[(h, c, int(s))] = lab
pieces = defaultdict(list)
for h, paf in (("hap1", paf1), ("hap2", paf2)):
    for l in open(paf):
        f = l.split("\t")
        if f[0] == "mid" or int(f[11]) < 20 or int(f[2]) > 100 or int(f[3]) < int(f[1]) - 100:
            continue
        a, b = int(f[7]) + FL, int(f[8]) - FL
        for t in seqs:
            lab, hh, rng = t.split("|")
            c, se = rng.rsplit(":", 1)
            s = int(se.split("-")[0]) - 1
            if hh == h and c == f[5] and s <= a and b <= s + len(seqs[t]):
                pieces[t].append((a - s, b - s, f[0].split("__")[0]))
rows = []
for t, s in seqs.items():
    for b in range(len(s) // BIN):
        v = cnt[t].get(b, [])
        if len(v) < BIN // 4:
            continue
        x = s[b * BIN:(b + 1) * BIN]
        acgt = sum(x.count(y) for y in "ACGT")
        gc = round(100 * (x.count("G") + x.count("C")) / max(1, acgt))
        rows.append((t, b, gc, statistics.median(v) / PEAK, het[(t, b)] / (BIN / 1000)))
ref = defaultdict(list)
for t, b, gc, cn, h in rows:
    if h >= 1.0:
        ref[gc].append(cn)
exp = {}
for g in range(20, 80):
    near = [v for d in (0, -1, 1, -2, 2) for v in ref.get(g + d, [])]
    if len(near) >= 8:
        exp[g] = statistics.median(near)
print("two-copy expectation (median k-mer CN of het bins) by GC %: " +
      " ".join(f"{g}:{exp[g]:.2f}" for g in sorted(exp) if g % 2 == 0))
print("window  bin_kb  GC  CN  het/kb  obs/exp   [* = inside a piece]")
summ = defaultdict(list)
for t, b, gc, cn, h in rows:
    if gc not in exp:
        continue
    r = cn / exp[gc]
    inp = [p for p in pieces.get(t, []) if p[0] <= b * BIN and (b + 1) * BIN <= p[1]]
    key = (t.split("|")[0], "piece" if inp else ("hetfree" if h < 0.3 else "het"))
    summ[key].append(r)
for (w, k), v in sorted(summ.items()):
    print(f"{w:12s} {k:8s} bins {len(v):3d}  obs/exp median {statistics.median(v):.2f}  "
          f"(IQR {statistics.quantiles(v, n=4)[0]:.2f}-{statistics.quantiles(v, n=4)[2]:.2f})" if len(v) >= 4 else
          f"{w:12s} {k:8s} bins {len(v):3d}  obs/exp " + " ".join(f"{x:.2f}" for x in v))
