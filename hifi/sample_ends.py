#!/usr/bin/env python3
"""Build classify_ends.py input for a sample of gaps: left flank and right flank from the gapped assembly, piece = the fill
another assembly holds (gapfill.py coordinates), plus fishing words next to each gap (unique in the gapped assembly).
Picks N_ONE gaps called one_copy and N_TWO called two_copies by fill_cn.py for one fill source.
Writes OUT/fills.fa (>gap_call middle=A-B), OUT/anchors.txt.
SOURCE_FILLS.fa = gapfill.py fills.fa of the same run (piece = fill without its 1 kb contig ends).
Usage: sample_ends.py GAPPED.fa GAPS.tsv FILL_CN.tsv SOURCE SOURCE_FILLS.fa OUT_DIR [N_ONE] [N_TWO]"""
import os
import random
import subprocess
import sys

gfa, gtsv, cn, src, sfa, out = sys.argv[1:7]
N_ONE = int(sys.argv[7]) if len(sys.argv) > 7 else 25
N_TWO = int(sys.argv[8]) if len(sys.argv) > 8 else 15
FL, NEAR, K, STEP = 10500, 3000, 25, 20
os.makedirs(out, exist_ok=True)
comp = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(comp)[::-1]
calls = {}
for l in list(open(cn))[1:]:
    f = l.rstrip("\n").split("\t")
    if f[1] == src and (f[0], f[1]) not in calls:
        calls[(f[0], f[1])] = f[6]
rows = [l.rstrip("\n").split("\t") for l in open(gtsv)]
h, rows = rows[0], rows[1:]
ix = {k: i for i, k in enumerate(h)}
random.seed(5)
one = [r for r in rows if calls.get((r[0], src)) == "one_copy" and r[ix[f"{src}_status"]] == "spanned"]
two = [r for r in rows if calls.get((r[0], src)) == "two_copies" and r[ix[f"{src}_status"]] == "spanned"]
pick = [(r, "one") for r in random.sample(one, min(N_ONE, len(one)))] + \
       [(r, "two") for r in random.sample(two, min(N_TWO, len(two)))]


def fetch(fa, regs):
    o = subprocess.run(["samtools", "faidx", fa] + regs, capture_output=True, text=True, check=True).stdout
    seqs, buf = [], []
    for l in o.splitlines():
        if l.startswith(">"):
            if buf or seqs:
                seqs.append("".join(buf))
            buf = []
        else:
            buf.append(l)
    seqs.append("".join(buf))
    return [s.upper() for s in seqs]


fills, name = {}, None
for l in open(sfa):
    if l.startswith(">"):
        f = l[1:].strip().split("|")
        name = (f[0], f[1]); fills[name] = ""
    elif name:
        fills[name] += l.strip().upper()
cand = set()
with open(os.path.join(out, "fills.fa"), "w") as fh:
    for r, call in pick:
        c, a, b = r[1], int(r[2]), int(r[3])               # N run [a, b) 0-based
        L, R = fetch(gfa, [f"{c}:{max(1, a - FL + 1)}-{a}", f"{c}:{b + 1}-{b + FL}"])
        if "N" in L[-NEAR:] or "N" in R[:NEAR] or (r[0], src) not in fills:
            continue
        P = fills[(r[0], src)][1000:-1000]                 # gapfill.py fills are in the gapped assembly's orientation
        name = f"{r[0]}_{call}"
        fh.write(f">{name} middle={len(L)}-{len(L) + len(P)}\n{L + P + R}\n")
        for s in (L[-NEAR:], R[:NEAR]):
            for i in range(0, len(s) - K, STEP):
                k = s[i:i + K]
                if "N" not in k:
                    cand.add(min(k, rc(k)))
cf = os.path.join(out, "cand.fa")
with open(cf, "w") as fh:
    fh.write("".join(f">{i}\n{k}\n" for i, k in enumerate(sorted(cand))))
subprocess.run(["jellyfish", "count", "-m", str(K), "-C", "-s", str(4 * len(cand)), "-t", "16", "--if", cf,
                "-o", os.path.join(out, "cand.jf"), gfa], check=True)
d = subprocess.run(["jellyfish", "dump", "-c", os.path.join(out, "cand.jf")], capture_output=True, text=True, check=True).stdout
keep = [k for k, n in (x.split() for x in d.splitlines()) if n == "1"]
open(os.path.join(out, "anchors.txt"), "w").write("\n".join(keep) + "\n")
os.remove(os.path.join(out, "cand.jf"))
print(f"{len(pick)} gaps sampled; {len(keep)}/{len(cand)} unique fishing words")
