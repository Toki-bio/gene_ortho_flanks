#!/usr/bin/env python3
"""Genome-wide test of hifiasm's two-copy rule on its own output.

For every p_utg unitig >= MIN_LEN: length, GC, rd (hifiasm's depth tag), and where its reads ended in the final Hi-C
haplotype contigs (A lines): in BOTH haplotypes, in one only, or neither. If hifiasm keeps a unitig in both haplotypes
only when its reads reach ~0.85 x the homozygous coverage (log: primary contig coverage range [21, infinity]), the share
placed in both should switch on sharply near rd 21, and GC should matter only through rd.
Output: table by rd bin (and by GC class within rd bins) of length-weighted shares: both / one haplotype / neither.
Usage: utg_rule_test.py P_UTG.gfa P_UTG.noseq.gfa HAP1.noseq.gfa HAP2.noseq.gfa [MIN_LEN]"""
import sys
from collections import defaultdict

putg, putg_ns, h1, h2 = sys.argv[1:5]
MIN_LEN = int(sys.argv[5]) if len(sys.argv) > 5 else 20000
info = {}
for l in open(putg):
    if l[0] != "S":
        continue
    f = l.rstrip("\n").split("\t")
    s = f[2]
    if len(s) < MIN_LEN:
        continue
    rd = next((int(x[5:]) for x in f[3:] if x.startswith("rd:i:")), None)
    gc = 100 * (s.count("G") + s.count("C") + s.count("g") + s.count("c")) / len(s)
    info[f[1]] = (len(s), gc, rd)
reads = defaultdict(list)
for l in open(putg_ns):
    if l[0] == "A":
        f = l.split("\t", 6)
        if f[1] in info:
            reads[f[1]].append(f[4])
want = {r for v in reads.values() for r in v}


def hap_reads(p):
    s = set()
    for l in open(p):
        if l[0] == "A":
            r = l.split("\t", 6)[4]
            if r in want:
                s.add(r)
    return s


H1, H2 = hap_reads(h1), hap_reads(h2)
tab = defaultdict(lambda: [0, 0, 0, 0])          # (rdbin, gcbin) -> length in both, one, neither, n unitigs
for u, (L, gc, rd) in info.items():
    rs = reads[u]
    if not rs or rd is None:
        continue
    b = sum(1 for r in rs if r in H1 and r in H2) / len(rs)
    o = sum(1 for r in rs if (r in H1) != (r in H2)) / len(rs)
    cls = 0 if b >= 0.5 else (1 if o >= 0.5 else 2)
    rb = min(rd, 40)
    gb = "<42" if gc < 42 else ("42-48" if gc < 48 else ">=48")
    for key in ((rb, "all"), (rb, gb)):
        tab[key][cls] += L
        tab[key][3] += 1
print("rd\tGC\tunitigs\tMb\tshare_both\tshare_one_hap\tshare_neither")
for (rb, gb) in sorted(tab, key=lambda k: (k[0], k[1] != "all", k[1])):
    a, o, n, k = tab[(rb, gb)]
    T = a + o + n
    if k >= 5:
        print(f"{rb}{'+' if rb == 40 else ''}\t{gb}\t{k}\t{T / 1e6:.1f}\t{a / T:.2f}\t{o / T:.2f}\t{n / T:.2f}")
