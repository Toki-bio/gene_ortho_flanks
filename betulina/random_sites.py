#!/usr/bin/env python3
"""Random control sites on a gapped assembly, in the gaps.tsv format of gapfill.py (gap chrom start end n_len class),
so gap_flank_cn.py prep extracts their 10 kb flanks the same way as for real gaps. Sites are >= AWAY bp from any
listed gap and from sequence ends, on sequences >= MINLEN; zero-length 'gaps' (start = end).
Usage: random_sites.py GENOME.fa.fai GAPS.tsv N SEED [AWAY=100000] [MINLEN=5000000] > random.tsv"""
import random
import sys

fai, gtsv, N, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
AWAY = int(sys.argv[5]) if len(sys.argv) > 5 else 100000
MINLEN = int(sys.argv[6]) if len(sys.argv) > 6 else 5000000
ln = {l.split()[0]: int(l.split()[1]) for l in open(fai) if int(l.split()[1]) >= MINLEN}
gaps = {}
for l in list(open(gtsv))[1:]:
    f = l.split("\t")
    gaps.setdefault(f[1], []).append((int(f[2]), int(f[3])))
rng = random.Random(seed)
tot = sum(ln.values())
out = []
while len(out) < N:
    x = rng.randrange(tot)
    for c, L in ln.items():
        if x < L:
            break
        x -= L
    if x < AWAY or x > L - AWAY or any(a - AWAY < x < b + AWAY for a, b in gaps.get(c, [])):
        continue
    out.append((c, x))
print("gap\tchrom\tstart\tend\tn_len\tclass")
for i, (c, x) in enumerate(sorted(out)):
    print(f"rnd{i + 1:04d}\t{c}\t{x}\t{x}\t0\trandom")
