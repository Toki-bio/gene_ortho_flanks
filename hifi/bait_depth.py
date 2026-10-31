#!/usr/bin/env python3
"""Total read depth (all identities) in 2 kb windows along each region's bait segment of the assembly (bait.fa, the
region +- its flanks in the deposited hap1), from the fished reads. Gap regions vs deletion/control regions of the
same read set: two copies of a segment give about twice the depth of one.
Usage: bait_depth.py BAIT.fa REGIONS.tsv READS.fq OUT_DIR [THREADS]"""
import os
import statistics
import subprocess
import sys
from collections import defaultdict

bait, regions, reads, out = sys.argv[1:5]
threads = sys.argv[5] if len(sys.argv) > 5 else "8"
W = 2000
os.makedirs(out, exist_ok=True)
paf = os.path.join(out, "bait_vs_reads.paf")
if not os.path.exists(paf):
    subprocess.run(f"minimap2 -c -x map-hifi --secondary=no -t {threads} {bait} {reads} > {paf} 2> /dev/null",
                   shell=True, check=True)
kind = {l.split("\t")[0]: l.split("\t")[1] for l in open(regions).read().splitlines()[1:]}
nrun = {}                                   # N-run spans per bait sequence
name, buf = None, []
for l in list(open(bait)) + [">"]:
    if l.startswith(">"):
        if name:
            s = "".join(buf)
            i, runs = 0, []
            while True:
                a = s.find("N", i)
                if a < 0:
                    break
                b = a
                while b < len(s) and s[b] == "N":
                    b += 1
                runs.append((a, b)); i = b
            nrun[name] = (len(s), runs)
        name, buf = l[1:].split()[0] if len(l) > 1 else None, []
    else:
        buf.append(l.strip().upper())
cov = defaultdict(lambda: defaultdict(set))
for l in open(paf):
    f = l.split("\t")
    t, ts, te = f[5], int(f[7]), int(f[8])
    if int(f[10]) < 2000:
        continue
    for w in range(-(-ts // W), te // W):
        cov[t][w].add(f[0])
for t in sorted(nrun):
    L, runs = nrun[t]
    reg = t.split("|")[0]
    depths = [len(cov[t][w]) for w in range(L // W)]
    gaps = ",".join(f"{a // 1000}-{b // 1000}kb" for a, b in runs) or "-"
    med = statistics.median(depths) if depths else 0
    print(f"{t}\t{kind.get(reg, '?')}\tlen {L // 1000} kb\tN-runs {gaps}\tmedian depth {med}\t" +
          " ".join(map(str, depths)))
