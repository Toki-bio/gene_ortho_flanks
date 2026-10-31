#!/usr/bin/env python3
"""GC-calibrated k-mer copy number of every gap fill found by wgs/gapfill.py: carried by both copies (fill it into the
gap) or by one copy (a real heterozygous insertion).

FILLS.fa: gapfill.py fills (>gap|hap|contig:lo-hi|strand|fill=N, fill with FLANK_OUT = 1 kb of contig on each side).
WIG: meryl-lookup -wig-count of FILLS.fa against the HiFi meryl db. CALIB: gc_calib.py output (its first line holds the
two-copy expectation per GC %). Copy number = median count of the fill's 21-mers (count > 0) / PEAK; obs/exp against
the expectation at the fill's GC. Fills shorter than MIN_LEN are not scored.
Usage: fill_cn.py FILLS.fa WIG CALIB [PEAK] > fill_cn.tsv"""
import statistics
import sys
from collections import defaultdict

fa, wig, calib = sys.argv[1:4]
PEAK = float(sys.argv[4]) if len(sys.argv) > 4 else 26
FO, MIN_LEN = 1000, 2000
exp = {}
for x in open(calib).readline().split(":", 1)[1].split():
    g, v = x.split(":")
    exp[int(g)] = float(v)
gs = sorted(exp)


def expected(gc):
    g = min(gs, key=lambda k: abs(k - gc))
    return exp[g]


seqs, name, buf = {}, None, []
for l in list(open(fa)) + [">"]:
    if l.startswith(">"):
        if name:
            seqs[name] = "".join(buf).upper()
        name, buf = l[1:].split()[0] if len(l) > 1 else None, []
    else:
        buf.append(l.strip())
cnt, chrom = defaultdict(list), None
for l in open(wig):
    if l.startswith("variableStep"):
        chrom = l.split("chrom=")[1].split()[0]; continue
    p, c = l.split()
    p, c = int(p), int(c)
    s = seqs.get(chrom)
    if s and FO <= p < len(s) - FO and c > 0:
        cnt[chrom].append(c)
print("gap\tsource\tfill_len\tGC\tcn\tobs_exp\tcall")
for n, s in seqs.items():
    gap, src = n.split("|")[:2]
    core = s[FO:len(s) - FO]
    if len(core) < MIN_LEN or len(cnt[n]) < len(core) // 4:
        print(f"{gap}\t{src}\t{len(core)}\t\t\t\tshort_or_unscored"); continue
    acgt = sum(core.count(x) for x in "ACGT")
    gc = 100 * (core.count("G") + core.count("C")) / max(1, acgt)
    cn = statistics.median(cnt[n]) / PEAK
    r = cn / expected(round(gc))
    call = "two_copies" if r >= 0.8 else ("one_copy" if r <= 0.65 else "unclear")
    if r >= 2.0:
        call = "repeat"
    print(f"{gap}\t{src}\t{len(core)}\t{gc:.1f}\t{cn:.2f}\t{r:.2f}\t{call}")
