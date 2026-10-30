#!/usr/bin/env python3
"""Why flanks were not placed: per flank, the first failing gate of gof.py.
Usage: why_unplaced.py OUT_DIR GENES.txt  (gof.py output dir; genes to report)"""
import collections
import os
import sys

out, genes_f = sys.argv[1:3]
genes = set(l.split()[0] for l in open(genes_f) if l.strip())
status = {}
for l in open(os.path.join(out, "verdicts.tsv")):
    f = l.rstrip("\n").split("\t")
    status[f[0]] = f[1]
flanks = {}
for l in open(os.path.join(out, "flanks.bed")):
    c, s, e, n = l.split()
    flanks[n] = (c, int(s), int(e))
pieces = collections.defaultdict(list)
for l in open(os.path.join(out, "flanks.paf")):
    f = l.split("\t")
    pieces[f[0]].append((int(f[3]) - int(f[2]), int(f[11]), "tp:A:P" in l, f[5]))

tab = collections.Counter()
contig_fail = collections.Counter()
for n, (c, s, e) in flanks.items():
    g, side, fl = n.split("|")
    if g not in genes:
        continue
    ql = e - s
    ps = [p for p in pieces.get(n, []) if p[2]]
    if not pieces.get(n):
        why = "no_hit"
    elif not any(p[1] >= 20 for p in ps):
        why = "only_mapq<20"
    else:
        cov = collections.Counter()
        for p in ps:
            if p[1] >= 20:
                cov[p[3]] += p[0]
        why = "placed" if max(cov.values()) >= min(0.5 * ql, 300) else "short_aligned"
    if ql < int(fl):
        why += "(flank_cut_by_contig_end)"
    tab[(status.get(g, "NA"), fl, why)] += 1
    if why != "placed":
        contig_fail[c] += 1
for k in sorted(tab):
    print(*k, tab[k], sep="\t")
print("== betulina contigs with most unplaced flanks")
for c, n in contig_fail.most_common(10):
    print(c, n)
