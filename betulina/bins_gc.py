#!/usr/bin/env python3
"""GC % per fixed bin of every sequence (streamed). Output: contig  start  end  GC%  ACGT_bp
Usage: bins_gc.py GENOME.fa BIN > gc.bed"""
import sys

fa, BIN = sys.argv[1], int(sys.argv[2])


def flush(name, s):
    for a in range(0, len(s), BIN):
        x = s[a:a + BIN]
        n = sum(x.count(c) for c in "ACGT")
        g = x.count("G") + x.count("C")
        print(f"{name}\t{a}\t{a + len(x)}\t{100 * g / n if n else -1:.2f}\t{n}")


name, buf = None, []
with open(fa) as fh:
    for l in fh:
        if l.startswith(">"):
            if name:
                flush(name, "".join(buf).upper())
            name, buf = l[1:].split()[0], []
        else:
            buf.append(l.strip())
if name:
    flush(name, "".join(buf).upper())
