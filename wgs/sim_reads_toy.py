#!/usr/bin/env python3
"""Simulated HiFi reads (15x, ~17 kb, random strand, 0.1 % substitutions) from FASTA contigs.
Usage: sim_reads_toy.py OUT.fq IN.fa [IN.fa ...]"""
import random
import sys

random.seed(9)
comp = str.maketrans("ACGT", "TGCA")
seqs = {}
for fn in sys.argv[2:]:
    name = None
    for l in open(fn):
        l = l.strip()
        if l.startswith(">"):
            name = fn + ":" + l[1:]
            seqs[name] = []
        else:
            seqs[name].append(l)
with open(sys.argv[1], "w") as h:
    k = 0
    for n, parts in seqs.items():
        s = "".join(parts)
        for _ in range(15 * len(s) // 17000):
            L = min(len(s), max(5000, int(random.gauss(17000, 3000))))
            p = random.randint(0, len(s) - L)
            r = list(s[p:p + L])
            for i in range(len(r)):
                if random.random() < 0.001:
                    r[i] = random.choice("ACGT")
            r = "".join(r)
            if random.random() < 0.5:
                r = r.translate(comp)[::-1]
            k += 1
            h.write(f"@r{k}_{n.split('/')[-1].replace('.fa:', '_')}\n{r}\n+\n{'I' * len(r)}\n")
