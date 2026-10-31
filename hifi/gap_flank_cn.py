#!/usr/bin/env python3
"""Per gap: GC-calibrated k-mer copy number of its two 10 kb flanks in the gapped assembly, next to the fill's copy
number (fill_cn.py) and the gapfill.py class. Two steps:
  prep:   gap_flank_cn.py prep GAPPED.fa GAPS.tsv OUT.fa          (flanks, stop at the next N run)
  table:  gap_flank_cn.py table OUT.fa WIG CALIB GAPS.tsv FILL_CN.tsv SOURCE [PEAK]
obs/exp ~1 = two copies, ~0.5 = one copy (expectation per GC from gc_calib.py's het bins)."""
import statistics
import subprocess
import sys
from collections import defaultdict

FL = 10000
if sys.argv[1] == "prep":
    gfa, gtsv, out = sys.argv[2:5]
    ln = {l.split()[0]: int(l.split()[1]) for l in open(gfa + ".fai")}
    regs, names = [], []
    for l in list(open(gtsv))[1:]:
        f = l.split("\t")
        c, a, b = f[1], int(f[2]), int(f[3])
        if a - FL >= 0:
            regs.append(f"{c}:{a - FL + 1}-{a}"); names.append(f"{f[0]}|L")
        if b + FL <= ln[c]:
            regs.append(f"{c}:{b + 1}-{b + FL}"); names.append(f"{f[0]}|R")
    o = subprocess.run(["samtools", "faidx", gfa] + regs, capture_output=True, text=True, check=True).stdout
    seqs, buf = [], None
    for l in o.splitlines():
        if l.startswith(">"):
            if buf is not None:
                seqs.append("".join(buf))
            buf = []
        else:
            buf.append(l)
    seqs.append("".join(buf))
    with open(out, "w") as fh:
        for n, s in zip(names, seqs):
            s = s.upper()
            # keep only the part next to the gap, up to the next N run
            s = s[s.rfind("N") + 1:] if n.endswith("|L") else s[:s.find("N")] if "N" in s else s
            if len(s) >= 3000:
                fh.write(f">{n}\n{s}\n")
    sys.exit(0)
fa, wig, calib, gtsv, cn, src = sys.argv[2:8]
PEAK = float(sys.argv[8]) if len(sys.argv) > 8 else 26
exp = {int(g): float(v) for g, v in (x.split(":") for x in open(calib).readline().split(":", 1)[1].split())}
expected = lambda gc: exp[min(exp, key=lambda k: abs(k - gc))]
seqs, name = {}, None
for l in open(fa):
    if l.startswith(">"):
        name = l[1:].strip(); seqs[name] = ""
    else:
        seqs[name] += l.strip()
cnt, cur = defaultdict(list), None
for l in open(wig):
    if l.startswith("variableStep"):
        cur = l.split("chrom=")[1].split()[0]; continue
    c = int(l.split()[1])
    if c > 0:
        cnt[cur].append(c)
fr = {}
for n, s in seqs.items():
    if len(cnt[n]) < len(s) // 4:
        continue
    gc = 100 * (s.count("G") + s.count("C")) / max(1, sum(s.count(x) for x in "ACGT"))
    fr[n] = statistics.median(cnt[n]) / PEAK / expected(round(gc))
fill = {}
for l in list(open(cn))[1:]:
    f = l.rstrip("\n").split("\t")
    if f[1] == src and f[0] not in fill:
        fill[f[0]] = (f[6], f[5])
cls = {l.split("\t")[0]: l.split("\t")[5] for l in list(open(gtsv))[1:]}
print("gap\tclass\tfill_call\tfill_obs_exp\tflankL_obs_exp\tflankR_obs_exp")
for g in cls:
    fc, fo = fill.get(g, ("-", ""))
    L, R = fr.get(f"{g}|L"), fr.get(f"{g}|R")
    print(f"{g}\t{cls[g]}\t{fc}\t{fo}\t{'' if L is None else f'{L:.2f}'}\t{'' if R is None else f'{R:.2f}'}")
