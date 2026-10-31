#!/usr/bin/env python3
"""K-mer copy number along windows: median HiFi count of the window's 21-mers per BIN (meryl-lookup -wig-count),
in units of the genome's two-copy peak (26 in ERR16751478, k=21). ~1.0 = two copies (homozygous or heterozygous: most
k-mers are shared), ~0.5 = one copy (hemizygous). Positions absent from the reads (count 0, e.g. assembly errors or
the other copy's SNPs) are reported as a fraction, not in the median.
Usage: kmer_cn.py WIG [PEAK] [BIN]"""
import statistics
import sys
from collections import defaultdict

wig = sys.argv[1]
PEAK = float(sys.argv[2]) if len(sys.argv) > 2 else 26
BIN = int(sys.argv[3]) if len(sys.argv) > 3 else 10000
cnt, chrom = defaultdict(lambda: defaultdict(list)), None
for l in open(wig):
    if l.startswith("variableStep"):
        chrom = l.split("chrom=")[1].split()[0]
        continue
    p, c = l.split()
    cnt[chrom][int(p) // BIN].append(int(c))
for t in cnt:
    bins = cnt[t]
    row = []
    for b in range(max(bins) + 1):
        v = bins.get(b, [])
        nz = [x for x in v if x > 0]
        row.append(f"{statistics.median(nz) / PEAK:.2f}" if len(nz) > BIN // 4 else "  - ")
    print(f"== {t}  copy number per {BIN // 1000} kb (two-copy peak = 1.00)")
    for i in range(0, len(row), 20):
        print(f"   {i * BIN // 1000:>5} kb  " + " ".join(row[i:i + 20]))
