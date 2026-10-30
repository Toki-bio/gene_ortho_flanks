#!/usr/bin/env python3
"""Why did the trizona assembly leave gaps where the reads rebuild sequence?
For each rebuilt gap span (spans/R0*.fa, header carries middle=a-b):
  A. middle aligned to the whole trizona assembly: bp found elsewhere (>= 95 % identity), and where
  B. read depth (all fished reads, map-hifi, MAPQ >= 0) in 1 kb windows: flanks vs middle; GC per window
  C. middle vs itself (off-diagonal) = internal repeats
Usage: whyfail.py QFA SPAN_DIR FISHED_FQ"""
import glob
import os
import statistics
import subprocess
import sys

qfa, sdir, fq = sys.argv[1:4]
WIN = 1000


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def read1(path):
    l = open(path).read().split("\n")
    a, b = map(int, l[0].split("middle=")[1].split("-"))
    return l[0][1:].split()[0], l[1], a, b


print("span\tmid_len\tA_elsewhere_bp\tA_top_hits\tB_depth_flank\tB_depth_mid_median\tB_depth_mid_min\tB_mid_low_windows\tB_gc_flank\tB_gc_mid\tB_gc_mid_max\tC_self_repeat_bp")
for f in sorted(glob.glob(os.path.join(sdir, "R0*.fa"))):
    name, seq, a, b = read1(f)
    mid = seq[a:b]
    d = os.path.dirname(f)
    mf = os.path.join(d, name + ".mid.fa")
    open(mf, "w").write(f">mid\n{mid}\n")
    # A: middle vs assembly
    cov, hits = set(), {}
    for l in run(["minimap2", "-c", "-x", "asm5", "-N", "20", "-t", "16", qfa, mf]).splitlines():
        x = l.split("\t")
        if int(x[9]) / int(x[10]) < 0.95:
            continue
        cov.update(range(int(x[2]), int(x[3])))
        k = x[5]
        hits[k] = hits.get(k, 0) + int(x[3]) - int(x[2])
    top = ",".join(f"{k}:{v}" for k, v in sorted(hits.items(), key=lambda t: -t[1])[:4])
    # B: depth along the span
    bam = os.path.join(d, name + ".bam")
    subprocess.run(f"minimap2 -ax map-hifi -t 16 {f} {fq} 2>/dev/null | samtools sort -@4 -o {bam} - && samtools index {bam}",
                   shell=True, check=True)
    dep = [0] * len(seq)
    for l in run(["samtools", "depth", "-a", "-Q", "0", bam]).splitlines():
        x = l.split("\t")
        dep[int(x[1]) - 1] = int(x[2])
    def wins(s, e):
        return [(statistics.mean(dep[i:i + WIN]),
                 sum(c in "GCgc" for c in seq[i:i + WIN]) / WIN) for i in range(s, e - WIN + 1, WIN)]
    fl = wins(0, a) + wins(b, len(seq))
    mi = wins(a, b)
    fd = statistics.median(w[0] for w in fl)
    md = [w[0] for w in mi]
    low = sum(1 for v in md if v < 0.5 * fd)
    # C: middle self-repeats (off-diagonal)
    rep = set()
    for l in run(["minimap2", "-c", "-x", "asm20", "-X", mf, mf]).splitlines():
        x = l.split("\t")
        if int(x[10]) >= 1000:
            rep.update(range(int(x[2]), int(x[3])))
    print("\t".join(map(str, [name, len(mid), len(cov), top or "-", f"{fd:.1f}", f"{statistics.median(md):.1f}",
          f"{min(md):.1f}", f"{low}/{len(md)}", f"{statistics.mean(w[1] for w in fl):.3f}",
          f"{statistics.mean(w[1] for w in mi):.3f}", f"{max(w[1] for w in mi):.3f}", len(rep)])))
