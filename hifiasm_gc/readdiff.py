#!/usr/bin/env python3
"""Reads gained/lost per haplotype (A lines) between asm25 and asm_gc, mapped to p_utg unitigs (rd, len)."""
import sys
from collections import Counter
A = "/home/toki/gof/asm25/tri.hic"
N = "/home/toki/gof/hifiasm_gc/asm_gc/tri.hic"
def rd(p):
    s = set()
    for l in open(p):
        if l[0] == "A":
            s.add(l.split("\t", 6)[4])
    return s
u_of = {}; info = {}
for l in open(A + ".p_utg.noseq.gfa"):
    if l[0] == "S":
        f = l.rstrip("\n").split("\t")
        r = next((x[5:] for x in f[3:] if x.startswith("rd:i:")), "?")
        L = next((x[5:] for x in f[3:] if x.startswith("LN:i:")), "?")
        info[f[1]] = (r, L)
    elif l[0] == "A":
        f = l.split("\t", 6); u_of[f[4]] = f[1]
for h in ("hap1", "hap2"):
    a = rd(f"{A}.{h}.p_ctg.noseq.gfa"); n = rd(f"{N}.{h}.p_ctg.noseq.gfa")
    for lab, s in (("gained", n - a), ("lost", a - n)):
        c = Counter(u_of.get(r, "-") for r in s)
        print(h, lab, len(s), "reads; top unitigs:", ", ".join(f"{u}(rd{info.get(u, ('?','?'))[0]},{info.get(u, ('?','?'))[1]}bp):{k}" for u, k in c.most_common(12)))
