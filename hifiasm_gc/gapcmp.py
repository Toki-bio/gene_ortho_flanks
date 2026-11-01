#!/usr/bin/env python3
"""Compare gapfill status of asm25 vs GC-corrected haplotypes per gap, by cause.
Usage: gapcmp.py GAPS.tsv CAUSES.tsv REPAIRS.tsv ASSEMBLY(dep1|dep2)"""
import sys
from collections import Counter, defaultdict
gaps, causes, repairs, asm = sys.argv[1:5]
cause = {}
for l in open(causes):
    f = l.rstrip("\n").split("\t")
    if f[0] != "gap":
        cause[f[0]] = f[13]
action = {}
for l in open(repairs):
    f = l.rstrip("\n").split("\t")
    if f[0] == asm:
        action[f[1]] = (f[5], f[6])
        cause.setdefault(f[1], f[5])
rows = []
hdr = None
for l in open(gaps):
    f = l.rstrip("\n").split("\t")
    if hdr is None:
        hdr = f; continue
    d = dict(zip(hdr, f))
    rows.append(d)
def sp(d, k):
    return d.get(k + "_status", "?") == "spanned"
tab = defaultdict(Counter)
cat1 = []
for d in rows:
    g = d["gap"]
    c = cause.get(g, "?")
    a = action.get(g, ("?", "?"))[1]
    grp = "CAT1_insert" if a == "insert" else c
    old = sp(d, "asm_h1") + sp(d, "asm_h2")
    new = sp(d, "gc_h1") + sp(d, "gc_h2")
    tab[grp][f"asm25:{old} gc:{new}"] += 1
    tab["ALL"][f"asm25:{old} gc:{new}"] += 1
    if a == "insert":
        cat1.append((g, d["chrom"], d["start"], old, new, d.get("asm_h1_status"), d.get("asm_h2_status"), d.get("gc_h1_status"), d.get("gc_h2_status"), d.get("gc_h1_fill", ""), d.get("gc_h2_fill", "")))
print("## haplotypes spanning the gap (0/1/2), asm25 vs GC-corrected, by group")
for grp in sorted(tab):
    print(grp, sum(tab[grp].values()), dict(sorted(tab[grp].items())), sep="\t")
print("## category-1 gaps")
print("gap\tchrom\tstart\tasm25_n\tgc_n\tasm_h1\tasm_h2\tgc_h1\tgc_h2\tgc_h1_fill\tgc_h2_fill")
for r in cat1:
    print("\t".join(map(str, r)))
