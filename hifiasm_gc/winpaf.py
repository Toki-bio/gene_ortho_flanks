#!/usr/bin/env python3
"""Per query window and per target assembly: best single-contig share of the window covered by alignments
(union of query intervals, primary+secondary, mapq any, >= 5 kb blocks). Usage: winpaf.py name=PAF [name=PAF ...]"""
import sys
from collections import defaultdict
res = defaultdict(dict)
qlen = {}
for arg in sys.argv[1:]:
    name, paf = arg.split("=", 1)
    iv = defaultdict(list)
    for l in open(paf):
        f = l.split("\t")
        if int(f[10]) < 5000:
            continue
        qlen[f[0]] = int(f[1])
        iv[(f[0], f[5])].append((int(f[2]), int(f[3])))
    best = defaultdict(lambda: (0, "-"))
    for (q, t), xs in iv.items():
        xs.sort(); cov = 0; s, e = xs[0]
        for a, b in xs[1:]:
            if a > e:
                cov += e - s; s, e = a, b
            else:
                e = max(e, b)
        cov += e - s
        if cov > best[q][0]:
            best[q] = (cov, t)
    for q in qlen:
        c, t = best[q]
        res[q][name] = (c / qlen[q], t)
names = [a.split("=", 1)[0] for a in sys.argv[1:]]
print("window\tlen\t" + "\t".join(f"{n}_best1contig\t{n}_contig" for n in names))
for q in sorted(qlen):
    print(q + "\t" + str(qlen[q]) + "\t" + "\t".join(f"{res[q].get(n, (0, '-'))[0]:.3f}\t{res[q].get(n, (0, '-'))[1]}" for n in names))
