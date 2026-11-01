#!/usr/bin/env python3
"""Contig ends of an assembly: 1 kb bins over the last SPAN bp of each end (for samtools bedcov), the last FL bp of each
end as FASTA (sequence features, mapping to other assemblies), and a 3 kb BED for read-overhang checks.
End names: contig|L (start) and contig|R (end). Bin name contig|L|k, k = kb from the contig end (0 = outermost).
A contig shorter than 2*SPAN gives each end at most half of it.
Usage: ends_prep.py GENOME.fa OUTPREFIX [SPAN=30000] [FL=20000] [MINCTG=0]"""
import subprocess
import sys

fa, out = sys.argv[1], sys.argv[2]
SPAN = int(sys.argv[3]) if len(sys.argv) > 3 else 30000
FL = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
MINCTG = int(sys.argv[5]) if len(sys.argv) > 5 else 0
ln = [(l.split()[0], int(l.split()[1])) for l in open(fa + ".fai")]
regs, names = [], []
with open(out + ".bins.bed", "w") as b, open(out + ".reads3k.bed", "w") as r:
    for c, L in ln:
        if L < MINCTG:
            continue
        half = min(SPAN, L // 2)
        for k in range(half // 1000):
            b.write(f"{c}\t{k * 1000}\t{(k + 1) * 1000}\t{c}|L|{k}\n")
            b.write(f"{c}\t{L - (k + 1) * 1000}\t{L - k * 1000}\t{c}|R|{k}\n")
        r.write(f"{c}\t0\t{min(3000, L)}\t{c}|L\n{c}\t{max(0, L - 3000)}\t{L}\t{c}|R\n")
        f = min(FL, L // 2)
        regs += [f"{c}:1-{f}", f"{c}:{L - f + 1}-{L}"]; names += [f"{c}|L", f"{c}|R"]
o = subprocess.run(["samtools", "faidx", fa] + regs, capture_output=True, text=True, check=True).stdout
seqs = [x.split("\n", 1)[1].replace("\n", "") for x in o.split(">")[1:]]
with open(out + ".ends.fa", "w") as fh:
    for n, s in zip(names, seqs):
        fh.write(f">{n}\n{s}\n")
