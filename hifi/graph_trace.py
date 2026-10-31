#!/usr/bin/env python3
"""Where did hifiasm put the reads of each gap region? (step 2 of the gap plan)

Reads of each region part (L flank, M piece, R flank) = reads with a >= 2 kb alignment to that part in
classify_ends.py's parts_vs_reads.paf (flank_profile.py writes it). For each graph (raw unitigs r_utg, cleaned
unitigs p_utg) prints the unitigs holding those reads with length, read depth (rd) and neighbours; then, for the
final haplotype contigs, how many of the part's reads ended in hap1 only / hap2 only / both / neither.
Usage: graph_trace.py PARTS_VS_READS.paf ASM_PREFIX(e.g. asm25/tri.hic)"""
import sys
from collections import Counter, defaultdict

paf, pre = sys.argv[1:3]
part_reads = defaultdict(set)
for l in open(paf):
    f = l.split("\t")
    if int(f[10]) >= 2000:
        part_reads[f[5]].add(f[0])
want = set().union(*part_reads.values())


def load(gfa):
    S, A, Lk = {}, {}, defaultdict(list)
    for l in open(gfa):
        f = l.rstrip("\n").split("\t")
        if f[0] == "S":
            S[f[1]] = (int(f[3][5:]), int(f[4][5:]) if len(f) > 4 and f[4].startswith("rd") else 0)
        elif f[0] == "A" and f[4] in want:
            A.setdefault(f[4], []).append(f[1])
        elif f[0] == "L":
            Lk[f[1]].append(f"{f[2]}>{f[3]}{f[4]}")
    return S, A, Lk


for g in ("r_utg", "p_utg"):
    S, A, Lk = load(f"{pre}.{g}.noseq.gfa")
    print(f"######## {g}: {len(S)} unitigs")
    for part in sorted(part_reads):
        rs = part_reads[part]
        c = Counter(u for r in rs for u in A.get(r, []))
        miss = sum(1 for r in rs if r not in A)
        top = ", ".join(f"{u}({n} reads; {S[u][0] // 1000} kb rd{S[u][1]})" for u, n in c.most_common(4))
        print(f"{part:24s} {len(rs):3d} reads, {miss} in no unitig | {top}")
        for u, n in c.most_common(2):
            nb = ", ".join(f"{x}" for x in Lk.get(u, [])[:6])
            print(f"{'':28s}{u} links: {nb}")
    if g == "p_utg":
        rd = sorted(v[1] for v in S.values() for _ in range(max(1, v[0] // 100000)))
        print(f"p_utg length-weighted median rd: {rd[len(rd) // 2]}")
hap = {}
for h in ("hap1", "hap2"):
    _, A, _ = load(f"{pre}.{h}.p_ctg.noseq.gfa")
    hap[h] = A
print("######## final haplotype contigs: part reads in hap1 only / hap2 only / both / neither  (top contigs)")
for part in sorted(part_reads):
    rs = part_reads[part]
    c = Counter(("h1" if r in hap["hap1"] else "") + ("h2" if r in hap["hap2"] else "") or "none" for r in rs)
    ctg = Counter(x for h in ("hap1", "hap2") for r in rs for x in hap[h].get(r, []))
    print(f"{part:24s} h1 {c['h1']:3d}  h2 {c['h2']:3d}  both {c['h1h2']:3d}  none {c['none']:3d} | "
          + ", ".join(f"{u}:{n}" for u, n in ctg.most_common(4)))
