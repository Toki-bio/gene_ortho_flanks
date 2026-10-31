#!/usr/bin/env python3
"""Sequences for characterising gaps: 10 kb flanks of the gaps no assembly spans (gapfill class split/unspanned), the
one-copy and two-copy fills (fill_cn.py calls, one source), and random control segments of the gapped assembly.
Also writes each gap's position along its sequence (distance to the nearest end, fraction of length) and a rodent
subset of the Dfam curated consensi for RepeatMasker -lib.
Usage: unknown_prep.py GAPPED.fa GAPS.tsv FILLS.fa FILL_CN.tsv SOURCE DFAM_CONSENSI DFAM_META OUT_DIR"""
import os
import random
import subprocess
import sys

gfa, gtsv, ffa, cn, src, dfa, dmeta, out = sys.argv[1:9]
FL, NCTRL = 10000, 60
os.makedirs(out, exist_ok=True)
ln = {l.split()[0]: int(l.split()[1]) for l in open(gfa + ".fai")}
rows = [l.rstrip("\n").split("\t") for l in open(gtsv)][1:]
unk = [r for r in rows if r[5] in ("split", "unspanned")]
regs, names = [], []
with open(os.path.join(out, "gap_positions.tsv"), "w") as fh:
    fh.write("gap\tchrom\tstart\tchrom_len\tdist_to_end\tfrac\tn_len\n")
    for r in unk:
        c, a, b = r[1], int(r[2]), int(r[3])
        fh.write(f"{r[0]}\t{c}\t{a}\t{ln[c]}\t{min(a, ln[c] - b)}\t{a / ln[c]:.3f}\t{b - a}\n")
        if a - FL >= 0:
            regs.append(f"{c}:{a - FL + 1}-{a}"); names.append(f"UNK_{r[0]}_L")
        if b + FL <= ln[c]:
            regs.append(f"{c}:{b + 1}-{b + FL}"); names.append(f"UNK_{r[0]}_R")
random.seed(9)
big = [c for c, n in ln.items() if n > 20_000_000]
for i in range(NCTRL):
    c = random.choice(big); s = random.randrange(0, ln[c] - FL)
    regs.append(f"{c}:{s + 1}-{s + FL}"); names.append(f"CTRL_{i}")
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
calls = {}
for l in list(open(cn))[1:]:
    f = l.rstrip("\n").split("\t")
    if f[1] == src and f[0] not in calls:
        calls[f[0]] = f[6]
with open(os.path.join(out, "feat.fa"), "w") as fh:
    for n, s in zip(names, seqs):
        if s.upper().count("N") < len(s) // 2:
            fh.write(f">{n}\n{s}\n")
    name = None
    for l in open(ffa):
        if l.startswith(">"):
            f = l[1:].split("|")
            name = f[0] if f[1] == src and calls.get(f[0]) in ("one_copy", "two_copies") else None
            if name:
                fh.write(f">{'ONE' if calls[name] == 'one_copy' else 'TWO'}_{name}\n")
        elif name:
            fh.write(l)
# rodent-relevant Dfam families
anc = {"Mammalia", "Theria", "Eutheria", "Boreoeutheria", "Euarchontoglires", "Glires", "Rodentia", "Amniota",
       "Tetrapoda", "Vertebrata", "Gnathostomata", "Euteleostomi", "Sarcopterygii", "root", "Eukaryota"}
keep = set()
for l in open(dmeta):
    f = l.rstrip("\n").split("\t")
    if f[2] in anc or "Rodentia" in f[3] or "Glires" in f[3] or f[2] in ("Muroidea", "Myomorpha", "Dipodidae", "Muridae"):
        keep.add(f[0])
with open(os.path.join(out, "rodent_lib.fa"), "w") as fh:
    take = False
    for l in open(dfa):
        if l.startswith(">"):
            take = l[1:].split()[0] in keep
        if take:
            fh.write(l)
print(f"{len(unk)} unspanned gaps; {len(keep)} Dfam families kept")
