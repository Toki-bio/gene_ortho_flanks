#!/usr/bin/env python3
"""Per contig end: reads reaching the end and reads OVERHANGING it (aligned to within SL bp of the end and soft/hard
clipped >= CLIP bp beyond it), plus where the overhanging part aligns (SA tag) - another contig's end means the
reads join two ends that the assembler left apart.
Input: SAM records (samtools view -q 20 -M -L reads3k.bed BAM) on stdin. Output per end:
end  ctg_len  n_reach  n_overhang  partner(best SA end, contig|L/R or contig:pos)  n_partner  partner_near_end(0/1)
Usage: samtools view ... | ends_reads.py GENOME.fa.fai"""
import re
import sys
from collections import Counter, defaultdict

SL, CLIP, NEAR = 500, 1000, 20000
ln = {l.split()[0]: int(l.split()[1]) for l in open(sys.argv[1])}
cig = re.compile(r"(\d+)([MIDNSHP=X])")
reach, over = Counter(), Counter()
part = defaultdict(Counter)
seen = set()
for l in sys.stdin:
    f = l.split("\t", 12)
    flag = int(f[1])
    if flag & 0x904 and not flag & 0x800:
        continue  # unmapped/secondary; supplementary kept (it is the piece of the read that is here)
    c, pos = f[2], int(f[3]) - 1
    ops = cig.findall(f[5])
    lc = int(ops[0][0]) if ops[0][1] in "SH" else 0
    rc = int(ops[-1][0]) if ops[-1][1] in "SH" else 0
    alen = sum(int(n) for n, o in ops if o in "MDN=X")
    end = pos + alen
    L = ln[c]
    sa = re.search(r"\tSA:Z:(\S+)", l)
    for side, near, clip in (("L", pos <= SL, lc), ("R", end >= L - SL, rc)):
        if not near:
            continue
        key = (f[0], c, side, flag & 0x800)
        if key in seen:
            continue
        seen.add(key)
        e = f"{c}|{side}"
        reach[e] += 1
        if clip >= CLIP:
            over[e] += 1
            if sa:
                for s in sa.group(1).rstrip(";").split(";"):
                    t, p = s.split(",")[:2]
                    if t == c and abs(int(p) - (pos if side == "L" else end)) < 50000:
                        continue
                    p, tl = int(p), ln[t]
                    nm = f"{t}|L" if p < NEAR else f"{t}|R" if p > tl - NEAR else f"{t}:{p // 10000 * 10}k"
                    part[e][nm] += 1
print("end\tctg_len\tn_reach\tn_overhang\tpartner\tn_partner\tpartner_near_end")
for c, L in ln.items():
    for side in "LR":
        e = f"{c}|{side}"
        if part[e]:
            p, n = part[e].most_common(1)[0]
        else:
            p, n = "-", 0
        print(f"{e}\t{L}\t{reach[e]}\t{over[e]}\t{p}\t{n}\t{int('|' in p)}")
