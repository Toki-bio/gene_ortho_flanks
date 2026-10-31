#!/usr/bin/env python3
"""Copy-number stretch and deletion junctions in windows of one haplotype (step 1b of the gap plan).

Reads (both copies) are aligned to WINDOWS.fa (minimap2 map-hifi, all chains). Prints per window:
  depth   reads with a >= 2 kb alignment covering each BIN (5 kb) bin (two copies ~ twice one copy)
  joins   read-supported jumps: two alignments of one read, adjacent in the read (<= 300 bp apart in the read), on the
          same window and strand, that skip >= MIN_JUMP bp of the window = the read's copy LACKS that stretch
          (clustered: same start and end +- 500 bp; reported with read counts)
Usage: sv_window.py WINDOWS.fa READS.fa|fq OUT_DIR [THREADS]"""
import os
import subprocess
import sys
from collections import defaultdict

win, reads, out = sys.argv[1:4]
threads = sys.argv[4] if len(sys.argv) > 4 else "16"
BIN, MIN_JUMP, SLACK = 5000, 2000, 300
os.makedirs(out, exist_ok=True)
paf = os.path.join(out, "win_vs_reads.paf")
if not os.path.exists(paf):
    subprocess.run(f"minimap2 -c -x map-hifi -t {threads} {win} {reads} > {paf} 2> /dev/null", shell=True, check=True)
CLIP = 2000
tlen, qlen, segs, depth = {}, {}, defaultdict(list), defaultdict(lambda: defaultdict(set))
for l in open(paf):
    f = l.split("\t")
    if "tp:A:S" in l:                     # secondary chains: repeats, not a second copy
        continue
    q, ql, qs, qe, st, t, tl, ts, te = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8])
    tlen[t] = tl
    qlen[q] = ql
    if te - ts >= 2000:
        for b in range(-(-ts // BIN), te // BIN):
            depth[t][b].add(q)
    fs, fe = (qs, qe) if st == "+" else (ql - qe, ql - qs)
    segs[(q, t)].append((fs, fe, ts, te, st))
joins = defaultdict(list)
for (q, t), sg in segs.items():
    sg.sort()
    for a, b in zip(sg, sg[1:]):
        if a[4] != b[4] or b[0] - a[1] > SLACK or a[1] - b[0] > SLACK:
            continue
        # in the window's frame: + strand a then b; - strand the read runs backwards along the window
        x, y = (a[3], b[2]) if a[4] == "+" else (b[3], a[2])
        if y - x >= MIN_JUMP:
            joins[t].append((x, y, q))
for t in sorted(tlen):
    d = [len(depth[t][b]) for b in range(tlen[t] // BIN)]
    print(f"== {t}  {tlen[t] // 1000} kb  depth per {BIN // 1000} kb: " + " ".join(map(str, d)))
    cl = []
    for x, y, q in sorted(joins[t]):
        for c in cl:
            if abs(c[0] - x) <= 500 and abs(c[1] - y) <= 500:
                c[2].add(q); break
        else:
            cl.append([x, y, {q}])
    for x, y, qs in sorted(cl, key=lambda c: -len(c[2])):
        if len(qs) >= 2:
            print(f"   join {x}-{y}  skips {y - x} bp  reads {len(qs)}")
    # clips: >= CLIP read bases beyond an alignment end that no other alignment of the read covers, away from the
    # window ends = the read's copy turns into sequence that is not in this window (divergent or elsewhere)
    cp = defaultdict(set)
    for (q, tt), sg in segs.items():
        if tt != t:
            continue
        ql = qlen[q]
        for s in sg:
            for side, rest, pos in (("R", ql - s[1], s[3] if s[4] == "+" else s[2]),
                                    ("L", s[0], s[2] if s[4] == "+" else s[3])):
                if rest < CLIP or pos < 5000 or pos > tlen[t] - 5000:
                    continue
                a, b = (s[1], s[1] + CLIP) if side == "R" else (s[0] - CLIP, s[0])
                if any(o is not s and o[0] < b and o[1] > a for o in sg):
                    continue
                # which way the window continues past the clip, in window coordinates
                wdir = ("+" if (side == "R") == (s[4] == "+") else "-")
                cp[(pos // 500, wdir)].add(q)
    for (b, wdir), qs in sorted(cp.items()):
        if len(qs) >= 3:
            print(f"   clip ~{b * 500}  window side {wdir}  reads {len(qs)}")
