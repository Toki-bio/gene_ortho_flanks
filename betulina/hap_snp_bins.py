#!/usr/bin/env python3
"""Het proxy for an assembly with two deposited haplotypes: SNPs between hap1 (target) and hap2 (query) per hap1 bin,
from minimap2 -x asm5 --cs PAF on stdin. Only primary alignments (tp:A:P) >= MIN_ALN bp. Output per bin:
contig  bin_start  aligned_bp  n_snp     (SNP/kb = n_snp / aligned_bp * 1000)
Usage: minimap2 -cx asm5 --cs hap1.fa hap2.fa | hap_snp_bins.py BIN > hapsnp.tsv"""
import re
import sys
from collections import defaultdict

BIN, MIN_ALN = int(sys.argv[1]), 50000
tok = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+)")
aln = defaultdict(int); snp = defaultdict(int)
for l in sys.stdin:
    f = l.rstrip("\n").split("\t")
    if "tp:A:P" not in f[12:] or int(f[8]) - int(f[7]) < MIN_ALN:
        continue
    t, p = f[5], int(f[7])
    cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
    for m in tok.findall(cs):
        if m[0] == ":":
            n = int(m[1:])
            while n > 0:  # aligned bp per bin
                b = p // BIN * BIN; k = min(n, b + BIN - p)
                aln[(t, b)] += k; p += k; n -= k
        elif m[0] == "*":
            b = p // BIN * BIN
            aln[(t, b)] += 1; snp[(t, b)] += 1; p += 1
        elif m[0] == "-":
            p += len(m) - 1
for (t, b), a in sorted(aln.items()):
    print(f"{t}\t{b}\t{a}\t{snp[(t, b)]}")
