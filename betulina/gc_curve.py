#!/usr/bin/env python3
"""Read depth vs GC, calibrated on two-copy bins (bins with heterozygous sites are certainly two copies).

Inputs (all per BIN, same bin starts): DEPTH.bed (contig start end depth; mosdepth --by or paf_depth.py),
GC.bed (bins_gc.py), HET (het_bins.py: contig bin n_het  -> mode 'clair'; hap_snp_bins.py: contig bin aligned n_snp ->
mode 'hapsnp'). EXCLUDE = file with contig names to skip (X / half-depth contigs), '-' for none.
Two-copy bins: het >= 1 per kb, >= 80 % ACGT, contig >= MINLEN. REF = median depth of all two-copy bins.
Line 1 (machine-readable, same format as gc_calib.py): 'two-copy expectation (median depth/REF of het bins) by GC %: g:v ...'
then REF, then a table: GC  n_twocopy  rel_twocopy  n_all  rel_all  rel_hetfree.
Usage: gc_curve.py DEPTH.bed GC.bed HET MODE EXCLUDE [MINLEN]"""
import statistics
import sys
from collections import defaultdict

dep, gcf, hetf, mode, excl = sys.argv[1:6]
MINLEN = int(sys.argv[6]) if len(sys.argv) > 6 else 1000000
ex = set() if excl == "-" else {l.split()[0] for l in open(excl) if l.strip()}
gc, clen = {}, defaultdict(int)
for l in open(gcf):
    c, s, e, g, n = l.split()
    clen[c] = max(clen[c], int(e))
    if float(g) >= 0 and int(n) >= 0.8 * (int(e) - int(s)) and int(e) - int(s) >= 1000:
        gc[(c, int(s))] = (round(float(g)), int(e) - int(s))
het = {}
for l in open(hetf):
    f = l.split()
    if mode == "clair":
        het[(f[0], int(f[1]))] = int(f[2])
    elif int(f[2]) > 0:
        het[(f[0], int(f[1]))] = (int(f[3]) / int(f[2]) * 1000, int(f[2]))
rows = []
for l in open(dep):
    c, s, e, d = l.split()[:4]
    k = (c, int(s))
    if c in ex or clen[c] < MINLEN or k not in gc:
        continue
    g, w = gc[k]
    if mode == "clair":
        h = het.get(k, 0) / w * 1000
    else:
        v = het.get(k)
        if v is None or v[1] < 0.8 * w:
            continue  # no 1:1 hap1-hap2 alignment over this bin: het state unknown
        h = v[0]
    rows.append((g, float(d), h))
two = [d for g, d, h in rows if h >= 1.0]
REF = statistics.median(two)
by = defaultdict(lambda: [[], [], []])
for g, d, h in rows:
    by[g][1].append(d / REF)
    if h >= 1.0:
        by[g][0].append(d / REF)
    elif h < 0.3:
        by[g][2].append(d / REF)
med = lambda v: f"{statistics.median(v):.2f}" if len(v) >= 30 else "NA"
exp = {}
for g in range(20, 80):
    near = [v for dd in (0, -1, 1, -2, 2) for v in by.get(g + dd, [[]])[0]]
    if len(near) >= 30:
        exp[g] = statistics.median(near)
print("two-copy expectation (median depth/REF of het bins) by GC %: " + " ".join(f"{g}:{exp[g]:.3f}" for g in sorted(exp)))
print(f"REF\t{REF:.2f}\tbins_used\t{len(rows)}\ttwo_copy_bins\t{len(two)}")
print("GC\tn_twocopy\trel_twocopy\tn_all\trel_all\trel_hetfree")
for g in sorted(by):
    t, a, f = by[g]
    if len(a) >= 30:
        print(f"{g}\t{len(t)}\t{med(t)}\t{len(a)}\t{med(a)}\t{med(f)}")
