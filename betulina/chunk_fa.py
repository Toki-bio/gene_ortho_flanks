#!/usr/bin/env python3
"""Split every sequence of a FASTA into SIZE bp chunks (named seq:start) so minimap2 can use all threads (it
parallelises over query sequences; one 200 Mb chromosome runs on one thread). Chunks that are all N are skipped.
Usage: chunk_fa.py IN.fa SIZE > chunks.fa"""
import sys

fa, SIZE = sys.argv[1], int(sys.argv[2])


def out(name, s):
    for a in range(0, len(s), SIZE):
        x = s[a:a + SIZE]
        if x.strip("Nn"):
            sys.stdout.write(f">{name}:{a}\n{x}\n")


name, buf = None, []
for l in open(fa):
    if l.startswith(">"):
        if name:
            out(name, "".join(buf))
        name, buf = l[1:].split()[0], []
    else:
        buf.append(l.strip())
if name:
    out(name, "".join(buf))
