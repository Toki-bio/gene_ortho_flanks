#!/usr/bin/env python3
"""Simulated HiFi reads with a known answer for rebuild.sh.

Truth per region (from regions.tsv):
  control       query as assembled                      -> rebuilt = assembly window
  gap, insert   N run replaced by the reference stretch -> rebuilt ~ stretch, r_cov ~ 1
  gap, empty    N run removed                           -> rebuilt ~ 0 bp
  deletion      query as assembled                      -> spanned, rebuilt = site, r_cov ~ 0
Reads: 20x, length N(17 kb, 4 kb) >= 3 kb, random strand, 0.1 % substitutions,
from region +- 60 kb; plus background reads from random query positions.

Usage: make_sim_toy.py regions.tsv QUERY_2BIT REF_2BIT OUT_DIR ID:truth [ID:truth ...]
  truth = control | insert | empty | asis
"""
import gzip
import os
import random
import subprocess
import sys

regions_f, q2bit, r2bit, out = sys.argv[1:5]
plan = dict(x.split(":") for x in sys.argv[5:])
FLANK, COV, BG = 60000, 20, 2000
random.seed(11)
os.makedirs(out, exist_ok=True)
comp = str.maketrans("ACGTN", "TGCAN")


def fetch(tb, c, s, e):
    r = subprocess.run(["twoBitToFa", f"{tb}:{c}:{s}-{e}", "stdout"], capture_output=True, text=True, check=True)
    return "".join(r.stdout.splitlines()[1:]).upper()


def mutate(s):
    s = list(s)
    for i in range(len(s)):
        if random.random() < 0.001:
            s[i] = random.choice("ACGT")
    return "".join(s)


reads = []
keep = []
for l in open(regions_f):
    f = l.rstrip("\n").split("\t")
    if f[0] == "id":
        keep.append(l)
        continue
    if f[0] not in plan:
        continue
    keep.append(l)
    qc, qs, qe, rc, rs, re_ = f[2], int(f[3]), int(f[4]), f[5], int(f[6]), int(f[7])
    left, mid, right = fetch(q2bit, qc, qs - FLANK, qs), fetch(q2bit, qc, qs, qe), fetch(q2bit, qc, qe, qe + FLANK)
    t = plan[f[0]]
    if t == "insert":
        mid = fetch(r2bit, rc, rs, re_)
    elif t == "empty":
        mid = ""
    seg = left + mid + right
    n = COV * len(seg) // 17000
    for i in range(n):
        L = max(3000, int(random.gauss(17000, 4000)))
        p = random.randint(0, max(0, len(seg) - L))
        r = mutate(seg[p:p + L])
        if random.random() < 0.5:
            r = r.translate(comp)[::-1]
        if "N" not in r:
            reads.append((f"{f[0]}_{t}_{i}", r))
    open(os.path.join(out, f"truth_{f[0]}.txt"), "w").write(f"{t}\t{len(mid)}\n")
# background reads from random query positions
chroms = [("OZ418341", 150_000_000), ("OZ418345", 150_000_000), ("OZ418350", 100_000_000)]
for i in range(BG):
    c, L0 = random.choice(chroms)
    p = random.randint(1_000_000, L0)
    r = fetch(q2bit, c, p, p + 17000)
    if "N" not in r:
        reads.append((f"bg_{i}", mutate(r)))
random.shuffle(reads)
with gzip.open(os.path.join(out, "reads.fastq.gz"), "wt") as h:
    for n, r in reads:
        h.write(f"@{n}\n{r}\n+\n{'I' * len(r)}\n")
open(os.path.join(out, "regions.tsv"), "w").write("".join(keep))
print(f"{len(reads)} reads -> {out}")
