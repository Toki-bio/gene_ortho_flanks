#!/usr/bin/env python3
"""Toy for orth_sites.py: flanks cut from the target itself, with known answers.
s1 point inside contig A (+)      -> spanned dist 0
s2 same, flanks reverse-complemented and swapped (site on - strand) -> spanned dist 0, strand -
s3 L = last 10 kb of contig A, R = first 10 kb of contig B -> split_end
s4 L from A, R from A 50 kb further (target 'insertion')    -> spanned dist 50000
s5 L from A, R = random sequence                            -> left_only
Usage: make_orth_toy.py A.fa B.fa > toy_flanks.fa  (writes toy_sites.tsv)"""
import random
import sys

def rd(p):
    return "".join(l.strip() for l in open(p) if not l.startswith(">")).upper()

A, B = rd(sys.argv[1]), rd(sys.argv[2])
rc = lambda s: s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
x = 2000000
fl = {"s1|L": A[x - 10000:x], "s1|R": A[x:x + 10000],
      "s2|L": rc(A[x + 1000000:x + 1010000]), "s2|R": rc(A[x + 990000:x + 1000000]),
      "s3|L": A[-10000:], "s3|R": B[:10000],
      "s4|L": A[x + 3000000 - 10000:x + 3000000], "s4|R": A[x + 3050000:x + 3060000],
      "s5|L": A[x + 4000000 - 10000:x + 4000000],
      "s5|R": "".join(random.Random(1).choice("ACGT") for _ in range(10000))}
for k, v in fl.items():
    print(f">{k}\n{v}")
with open("toy_sites.tsv", "w") as fh:
    fh.write("gap\n" + "".join(f"s{i}\n" for i in range(1, 6)))
