#!/usr/bin/env python3
"""Toy for gapfill.py with a known answer. From N-free 60 kb windows of the query genome:
the toy genome has segment X (length below) replaced by 100 N; the haplotype contigs carry
the window (0.1 % substitutions) with or without X. Truth per case (hap1, hap2 fill):
  insert   both carry X (20 kb)                  filled   20000 / 20000
  empty    X = 0 (100 N inserted, nothing lost)  joined   0 / 0
  minus    hap1 contig reverse-complemented      filled   8000 / 8000
  endsin   hap1 contig ends inside X, hap2 none  unspanned (left_only / none)
  split    hap1 window cut inside X into 2 contigs, hap2 none   split
  het      hap1 carries X (5 kb), hap2 lacks it  filled   5000 / 0
Usage: make_gapfill_toy.py QUERY_2BIT OUT_DIR"""
import os
import random
import subprocess
import sys

q2bit, out = sys.argv[1:3]
os.makedirs(out, exist_ok=True)
random.seed(5)
comp = str.maketrans("ACGT", "TGCA")


def fetch(c, s, e):
    r = subprocess.run(["twoBitToFa", f"{q2bit}:{c}:{s}-{e}", "stdout"], capture_output=True, text=True, check=True)
    return "".join(r.stdout.splitlines()[1:]).upper()


def mut(s):
    s = list(s)
    for i in range(len(s)):
        if random.random() < 0.001:
            s[i] = random.choice("ACGT")
    return "".join(s)


cases = [("insert", 20000), ("empty", 0), ("minus", 8000), ("endsin", 6000), ("split", 6000), ("het", 5000)]
genome, h1, h2, truth = [], [], [], []
pos = 30_000_000
for name, x in cases:
    while True:
        w = fetch("OZ418341", pos, pos + 60000)
        pos += 200000
        if "N" not in w:
            break
    m = 25000
    gapped = w[:m] + "N" * 100 + w[m + x:]
    genome.append((f"toy_{name}", gapped))
    full = mut(w)
    if name == "minus":
        h1.append((f"{name}_h1", full.translate(comp)[::-1]))
        h2.append((f"{name}_h2", mut(w)))
    elif name == "endsin":
        h1.append((f"{name}_h1", full[:m + x // 2]))
    elif name == "split":
        h1.append((f"{name}_h1a", full[:m + x // 2]))
        h1.append((f"{name}_h1b", full[m + x // 2:]))
    elif name == "het":
        h1.append((f"{name}_h1", full))
        h2.append((f"{name}_h2", mut(w[:m] + w[m + x:])))
    else:
        h1.append((f"{name}_h1", full))
        h2.append((f"{name}_h2", mut(w)))
    truth.append(name)


def write(fn, recs):
    with open(os.path.join(out, fn), "w") as h:
        for n, s in recs:
            h.write(f">{n}\n{s}\n")


write("genome.fa", genome)
write("hap1.fa", h1)
write("hap2.fa", h2)
print("toy written:", ", ".join(truth))
