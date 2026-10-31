#!/usr/bin/env python3
"""Toy oracle for polish_illumina.sh.
make TRUE.fa DAMAGED.fa SITES.tsv: damage TRUE.fa with homopolymer shortenings (runs >= 4 lose one base)
and substitutions at known positions.
score TRUE.fa POLISHED.fa: minimap2 -c --cs=long alignment, count remaining differences.
Usage: polish_toy.py make IN.fa TRUE.fa DAMAGED.fa SITES.tsv | polish_toy.py score TRUE.fa POLISHED.fa"""
import random
import re
import subprocess
import sys


def read1(fn):
    return "".join(l.strip() for l in open(fn) if not l.startswith(">")).upper()


def make(src, true_fa, dmg_fa, sites_fn, n_hp=200, n_sub=50):
    random.seed(5)
    s = read1(src)
    s = s[:2_000_000].replace("N", "")
    runs = [m.start() for m in re.finditer(r"A{4,}|C{4,}|G{4,}|T{4,}", s)]
    hp = sorted(random.sample(runs, n_hp))
    taken = set(hp)
    sub = []
    while len(sub) < n_sub:
        p = random.randrange(1000, len(s) - 1000)
        if all(abs(p - q) > 50 for q in taken):
            sub.append(p)
            taken.add(p)
    sub.sort()
    d = list(s)
    for p in sub:
        d[p] = random.choice([b for b in "ACGT" if b != s[p]])
    for p in hp:
        d[p] = ""
    with open(true_fa, "w") as h:
        h.write(">toy\n" + s + "\n")
    with open(dmg_fa, "w") as h:
        h.write(">toy\n" + "".join(d) + "\n")
    with open(sites_fn, "w") as h:
        for p in hp:
            h.write(f"hp_del\t{p}\n")
        for p in sub:
            h.write(f"sub\t{p}\n")
    print(f"true {len(s)} bp; damaged with {n_hp} homopolymer deletions + {n_sub} substitutions")


def score(true_fa, pol_fa):
    out = subprocess.run(["minimap2", "-c", "--cs=long", "-x", "asm5", true_fa, pol_fa],
                         capture_output=True, text=True, check=True).stdout
    best = max((l.split("\t") for l in out.splitlines()), key=lambda f: int(f[10]))
    cs = next(x[5:] for x in best[12:] if x.startswith("cs:Z:"))
    ops = re.findall(r"(=[A-Za-z]+|\*[a-z][a-z]|[-+][a-z]+)", cs)
    subs = sum(1 for o in ops if o[0] == "*")
    indel = sum(1 for o in ops if o[0] in "+-")
    print(f"aligned {best[7]}-{best[8]} of {best[6]} (true); remaining substitutions {subs}, indels {indel}")


if sys.argv[1] == "make":
    make(*sys.argv[2:6])
else:
    score(*sys.argv[2:4])
