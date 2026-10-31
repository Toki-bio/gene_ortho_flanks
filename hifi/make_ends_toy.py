#!/usr/bin/env python3
"""Toy truth for classify_ends.py. Haplotype A carries a piece between flanks L and R; haplotype B differs by case:
  hom_ins     B carries the same piece (0.3 % diverged)         -> both sides: piece only
  het_del     B joins L straight to R                           -> piece + direct
  het_other   B carries a different 5 kb sequence               -> piece + other
  het_drop    B has no reads within 8 kb of the gap             -> piece (+ ends)
  hom_join    A and B both join L to R (piece empty)            -> direct only
Reads: 15 kb, 0.1 % errors, 12x per haplotype, random strand. Writes OUT/fills.fa (>case middle=A-B) and OUT/reads.fq.
Usage: make_ends_toy.py OUT_DIR"""
import os
import random
import sys

random.seed(11)
out = sys.argv[1]
os.makedirs(out, exist_ok=True)
rs = lambda n: "".join(random.choice("ACGT") for _ in range(n))
rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def mut(s, r):
    return "".join(random.choice([c for c in "ACGT" if c != x]) if random.random() < r else x for x in s)


F = 30000
n = 0
with open(os.path.join(out, "fills.fa"), "w") as fa, open(os.path.join(out, "reads.fq"), "w") as fq:
    for case in ["hom_ins", "het_del", "het_other", "het_drop", "hom_join"]:
        L, R = rs(F), rs(F)
        P = "" if case == "hom_join" else rs(20000)
        A = L + P + R
        if case == "hom_ins":
            B = mut(L + P + R, 0.003)
        elif case == "het_other":
            B = mut(L, 0.003) + rs(5000) + mut(R, 0.003)
        else:
            B = mut(L, 0.003) + mut(R, 0.003)
        fa.write(f">{case} middle=10500-{10500 + len(P)}\n{L[-10500:] + P + R[:10500]}\n")
        for hap, s in (("A", A), ("B", B)):
            for _ in range(int(12 * len(s) / 15000)):
                st = random.randrange(0, len(s) - 15000)
                if case == "het_drop" and hap == "B" and st < F + 8000 and st + 15000 > F - 8000:
                    continue
                r = mut(s[st:st + 15000], 0.001)
                r = r if random.random() < 0.5 else rc(r)
                fq.write(f"@t{n}_{case}_{hap}\n{r}\n+\n{'I' * len(r)}\n")
                n += 1
print(f"{n} reads -> {out}")
