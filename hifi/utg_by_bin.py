#!/usr/bin/env python3
"""Per 50 kb bin of each window: the hifiasm unitigs (p_utg) its reads belong to, their rd, and where those reads ended
in the final haplotype contigs (hap1 / hap2 / both / neither). Calibrates hifiasm's rd and its homozygous handling on
the control windows (CTRL1 = run of homozygosity by het_window.py) against the gap regions.
Reads per bin: primary alignments >= 2 kb whose midpoint falls in the bin (win_vs_reads.paf of sv_window.py).
Usage: utg_by_bin.py WIN_VS_READS.paf ASM_PREFIX [BIN]"""
import sys
from collections import Counter, defaultdict

paf, pre = sys.argv[1:3]
BIN = int(sys.argv[3]) if len(sys.argv) > 3 else 50000
bins = defaultdict(set)
for l in open(paf):
    f = l.split("\t")
    if "tp:A:P" not in l or int(f[8]) - int(f[7]) < 2000:
        continue
    bins[(f[5], (int(f[7]) + int(f[8])) // 2 // BIN)].add(f[0])
want = set().union(*bins.values())


def load(gfa):
    S, A = {}, {}
    for l in open(gfa):
        if l[0] == "S":
            f = l.split("\t")
            S[f[1]] = (int(f[3][5:]), int(f[4][5:]))
        elif l[0] == "A":
            f = l.split("\t", 6)
            if f[4] in want:
                A.setdefault(f[4], set()).add(f[1])
    return S, A


S, U = load(f"{pre}.p_utg.noseq.gfa")
_, H1 = load(f"{pre}.hap1.p_ctg.noseq.gfa")
_, H2 = load(f"{pre}.hap2.p_ctg.noseq.gfa")
for t in sorted({t for t, _ in bins}):
    print(f"== {t}")
    for b in sorted(b for tt, b in bins if tt == t):
        rs = bins[(t, b)]
        c = Counter(u for r in rs for u in U.get(r, ()))
        top = " ".join(f"{u}:{n}r/{S[u][0] // 1000}kb/rd{S[u][1]}" for u, n in c.most_common(2))
        h = Counter(("1" if r in H1 else "") + ("2" if r in H2 else "") or "0" for r in rs if r in U)
        print(f"   {b * BIN // 1000:>4} kb  reads {len(rs):3d}  in utg {sum(1 for r in rs if r in U):3d} | hap1 {h['1']:3d} "
              f"hap2 {h['2']:3d} both {h['12']:3d} none {h['0']:3d} | {top}")
