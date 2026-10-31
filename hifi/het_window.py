#!/usr/bin/env python3
"""Heterozygous-site density, depth and GC per 50 kb bin of each window (tells one copy / two identical copies /
two different copies apart).
Reads are aligned with minimap2 -c --cs (primary chains >= 2 kb). A position is a het site when >= MIN_DEPTH reads
cover it, the most common substitution there is seen in >= MIN_ALT reads and makes up 25-75 % of the covering reads.
Expectations: two copies, normal heterozygosity -> ~3 het/kb (trizona hap1 vs hap2: 3.03 SNVs/kb); one copy
(hemizygous) or two identical copies (ROH) -> ~0; depth separates those two (one copy ~ half).
Usage: het_window.py WINDOWS.fa READS.fa|fq OUT_DIR [THREADS] [BIN]"""
import os
import re
import subprocess
import sys
from collections import defaultdict

win, reads, out = sys.argv[1:4]
threads = sys.argv[4] if len(sys.argv) > 4 else "16"
BIN = int(sys.argv[5]) if len(sys.argv) > 5 else 50000
MIN_DEPTH, MIN_ALT = 6, 3
os.makedirs(out, exist_ok=True)
paf = os.path.join(out, "win_vs_reads.cs.paf")
if not os.path.exists(paf):
    subprocess.run(f"minimap2 -c --cs -x map-hifi --secondary=no -t {threads} {win} {reads} > {paf} 2> /dev/null",
                   shell=True, check=True)
seqs, name, buf = {}, None, []
for l in list(open(win)) + [">"]:
    if l.startswith(">"):
        if name:
            seqs[name] = "".join(buf).upper()
        name, buf = l[1:].split()[0] if len(l) > 1 else None, []
    else:
        buf.append(l.strip())
cov = {t: [0] * (len(s) + 1) for t, s in seqs.items()}         # difference array
subs = defaultdict(lambda: defaultdict(int))                   # (t, pos) -> alt base -> count
tok = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+)")
for l in open(paf):
    f = l.rstrip("\n").split("\t")
    if int(f[8]) - int(f[7]) < 2000:
        continue
    t, ts, te = f[5], int(f[7]), int(f[8])
    cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
    cov[t][ts] += 1; cov[t][te] -= 1
    p = ts
    for m in tok.findall(cs):
        if m[0] == ":":
            p += int(m[1:])
        elif m[0] == "*":
            subs[(t, p)][m[2].upper()] += 1
            p += 1
        elif m[0] == "-":
            p += len(m) - 1
for t in sorted(seqs):
    c, run = [], 0
    for x in cov[t][:-1]:
        run += x; c.append(run)
    s = seqs[t]
    print(f"== {t}  {len(s) // 1000} kb   per {BIN // 1000} kb bin: het/kb  depth  GC%")
    rows = []
    for b in range(len(s) // BIN):
        a, z = b * BIN, (b + 1) * BIN
        het = 0
        for p in range(a, z):
            d = subs.get((t, p))
            if d and c[p] >= MIN_DEPTH:
                n = max(d.values())
                if n >= MIN_ALT and 0.25 <= n / c[p] <= 0.75:
                    het += 1
        dep = sorted(c[a:z])[BIN // 2]
        acgt = sum(s[a:z].count(x) for x in "ACGT")
        gc = 100 * sum(s[a:z].count(x) for x in "GC") / max(1, acgt)
        rows.append(f"{het / (BIN / 1000):.2f}/{dep}/{gc:.0f}")
    print("   " + "  ".join(rows))
