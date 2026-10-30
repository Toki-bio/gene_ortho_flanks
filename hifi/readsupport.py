#!/usr/bin/env python3
"""Reads per haplotype contig: all spanning contigs of a region (spans/<rid>__*.fa) are one
target, so each read goes to the haplotype it matches best; a read supports a contig when
its primary alignment covers the whole middle plus 500 bp each side (left/right = covers
that junction only). Also writes, per contig, a SAM (qualities dropped) for ViewAlign.
Usage: readsupport.py SPAN_DIR ASM_DIR OUT_DIR RID [RID ...]"""
import glob
import os
import subprocess
import sys

sdir, adir, odir = sys.argv[1:4]
M = 500
os.makedirs(odir, exist_ok=True)
print("region\tcontig\tmiddle_len\treads_spanning\treads_left_only\treads_right_only\treads_total")
for rid in sys.argv[4:]:
    fas = sorted(f for f in glob.glob(os.path.join(sdir, f"{rid}__*.fa")) if not f.endswith(".mid.fa"))
    mids, ref = {}, os.path.join(odir, f"{rid}.ref.fa")
    with open(ref, "w") as h:
        for f in fas:
            l = open(f).read().split("\n")
            name = l[0][1:].split()[0]
            a, b = map(int, l[0].split("middle=")[1].split("-"))
            mids[name] = (a, b)
            h.write(f">{name}\n{l[1]}\n")
            open(os.path.join(odir, f"{name}.fa"), "w").write(f">{name}\n{l[1]}\n")
    sam = subprocess.run(f"minimap2 -ax map-hifi -t 16 --secondary=no {ref} {adir}/{rid}/reads.fq 2>/dev/null",
                         shell=True, capture_output=True, text=True, check=True).stdout
    cnt = {n: [0, 0, 0, 0] for n in mids}
    per = {n: [] for n in mids}
    hdr = [l for l in sam.splitlines() if l.startswith("@HD")]
    for l in sam.splitlines():
        if l.startswith("@"):
            continue
        x = l.split("\t")
        flag = int(x[1])
        if flag & 0x904 or x[2] not in mids:
            continue
        s = int(x[3]) - 1
        span = 0
        num = ""
        for ch in x[5]:
            if ch.isdigit():
                num += ch
            else:
                if ch in "MDN=X":
                    span += int(num)
                num = ""
        e = s + span
        a, b = mids[x[2]]
        L, R = s <= a - M and e >= a + min(M, (b - a) // 2 or 1), s <= b - min(M, (b - a) // 2 or 1) and e >= b + M
        c = cnt[x[2]]
        c[3] += 1
        if s <= a - M and e >= b + M:
            c[0] += 1
        elif L:
            c[1] += 1
        elif R:
            c[2] += 1
        x[10] = "*"
        per[x[2]].append("\t".join(x[:11]))
    for n, (a, b) in mids.items():
        c = cnt[n]
        print(f"{rid}\t{n}\t{b - a}\t{c[0]}\t{c[1]}\t{c[2]}\t{c[3]}")
        seqlen = len(open(os.path.join(odir, f"{n}.fa")).read().split("\n")[1])
        with open(os.path.join(odir, f"{n}.sam"), "w") as h:
            h.write(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{n}\tLN:{seqlen}\n")
            h.write("\n".join(sorted(per[n], key=lambda r: int(r.split("\t")[3]))) + "\n")
