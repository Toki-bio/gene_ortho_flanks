#!/usr/bin/env python3
"""Which haplotype does each read near the gap belong to? Each read is mapped to hap1 and to hap2 alone; its mismatch rate on the
two flanks (hap1 and hap2 differ there) assigns it. Then: do hap1-assigned reads cross the gap, and into what?"""
import subprocess, re, sys, collections
S = "/home/toki/miniforge3/envs/sinederella/bin/samtools"
G = "/home/toki/gof/gene_views"


def load(p):
    R = {}
    for l in subprocess.run([S, "view", "-F", "0x800", p], capture_output=True, text=True).stdout.splitlines():
        f = l.split("\t")
        nm = int(re.search(r"NM:i:(\d+)", l).group(1))
        ops = re.findall(r"(\d+)([MIDNSH=X])", f[5])
        pos = int(f[3]) - 1
        ref = sum(int(a) for a, b in ops if b in "MDN=X")
        aln = sum(int(a) for a, b in ops if b in "M=X")
        big = sum(int(a) for a, b in ops if b in "ID" and int(a) >= 50)
        small = nm - big
        clipL = int(ops[0][0]) if ops[0][1] in "SH" else 0
        clipR = int(ops[-1][0]) if ops[-1][1] in "SH" else 0
        bigI = max([int(a) for a, b in ops if b == "I"] + [0])
        R[f[0]] = dict(s=pos, e=pos + ref, aln=aln, rate=small / max(aln, 1), nm=small, cl=clipL, cr=clipR, bigI=bigI, L=len(f[9]))
    return R


for g, g0, g1 in (("TMED4", 25000, 25100), ("EPS8L1", 25000, 25100)):
    h1 = load(f"{G}/{g}/{g}_tri_hap1.bam"); h2 = load(f"{G}/{g}/{g}_tri_hap2.bam")
    print("==", g, "reads on hap1", len(h1), "on hap2", len(h2))
    near = {k: v for k, v in h1.items() if v["s"] < g1 + 200 and v["e"] > g0 - 200}
    c = collections.Counter()
    for k, v in sorted(near.items(), key=lambda x: x[1]["s"]):
        w = h2.get(k)
        a = v["nm"]; b = w["nm"] if w else None
        call = "hap?" if w is None else "hap1" if a + 3 < b else "hap2" if b + 3 < a else "tie"
        cross = v["s"] <= g0 - 500 and v["e"] >= g1 + 500
        c[(call, "crosses" if cross else "stops")] += 1
        print(f" {k[-22:]} len {v['L']} hap1 {v['s']}-{v['e']} nm {a} clipL {v['cl']} clipR {v['cr']} bigI {v['bigI']} | hap2 {(w['s'], w['e'], b, w['cl'], w['cr']) if w else None} -> {call} {'CROSS' if cross else ''}")
    print(" summary", dict(c))
