#!/usr/bin/env python3
"""Per-exon presence in each row of an aln_rebuild.py alignment.
For every BED12 transcript on the reference row, and every exon, report for each other row:
aligned (non-gap) bases over the exon columns, and identity to the reference over those bases.
Usage: exon_presence.py ALN.fa ALN.bed"""
import sys

aln_f, bed_f = sys.argv[1:3]
rows, name = {}, None
for l in open(aln_f):
    l = l.rstrip()
    if l.startswith(">"):
        name = l[1:].split()[0]
        rows[name] = []
    else:
        rows[name].append(l.upper())
rows = {k: "".join(v) for k, v in rows.items()}


def col_of(seq):
    """ungapped position -> alignment column"""
    return [i for i, c in enumerate(seq) if c != "-"]


print("transcript\texon\tref_start\tref_len\trow\taligned\tident")
for l in open(bed_f):
    f = l.rstrip("\n").split("\t")
    if len(f) < 12 or f[0] not in rows:
        continue
    ref = rows[f[0]]
    cols = col_of(ref)
    start = int(f[1])
    sizes = [int(x) for x in f[10].strip(",").split(",")]
    offs = [int(x) for x in f[11].strip(",").split(",")]
    for k, (sz, of) in enumerate(zip(sizes, offs), 1):
        a, b = start + of, start + of + sz
        ec = cols[a:b]
        for r, s in rows.items():
            if r == f[0]:
                continue
            al = [(ref[c], s[c]) for c in ec if s[c] != "-"]
            idt = sum(x == y for x, y in al) / len(al) if al else 0
            print(f"{f[3]}\t{k}\t{a}\t{sz}\t{r}\t{len(al)}\t{idt:.3f}")
