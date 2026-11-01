#!/usr/bin/env python3
"""Heterozygous SNPs per bin from a Clair3 VCF (PASS, GT 0/1, SNV, AF 0.35-0.65 = the ~1 % false-het filter measured
on the hemizygous X). Output: contig  bin_start  n_het
Usage: zcat calls.vcf.gz | het_bins.py BIN > het.tsv"""
import sys
from collections import Counter

BIN = int(sys.argv[1])
c = Counter()
for l in sys.stdin:
    if l[0] == "#":
        continue
    f = l.split("\t")
    if f[6] != "PASS" or len(f[3]) != 1 or len(f[4]) != 1:
        continue
    fmt, smp = f[8].split(":"), f[9].rstrip("\n").split(":")
    d = dict(zip(fmt, smp))
    if d.get("GT") not in ("0/1", "0|1", "1|0"):
        continue
    try:
        af = float(d.get("AF", "nan"))
    except ValueError:
        continue
    if 0.35 <= af <= 0.65:
        p = int(f[1]) - 1
        c[(f[0], p // BIN * BIN)] += 1
for (k, b), n in sorted(c.items()):
    print(f"{k}\t{b}\t{n}")
