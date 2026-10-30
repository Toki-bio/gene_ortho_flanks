#!/usr/bin/env python3
"""GC and lowercase (soft-masked) fraction of betulina stretches that face a trizona
assembly gap, vs random betulina windows of the same size.
Usage: gap_regions_gc.py PROFILE_N.tsv REF_2BIT"""
import random, re, subprocess, sys
prof, r2bit = sys.argv[1:3]
def fetch(c, s, e):
    return "".join(subprocess.run(["twoBitToFa", f"{r2bit}:{c}:{s}-{e}", "stdout"], capture_output=True, text=True, check=True).stdout.splitlines()[1:])
def stats(seq):
    acgt = sum(seq.upper().count(b) for b in "ACGT") or 1
    gc = (seq.upper().count("G") + seq.upper().count("C")) / acgt
    low = sum(1 for ch in seq if ch.islower()) / max(1, len(seq))
    # simple repeats: share of bases in di/tri-nucleotide runs >= 20 bp
    simple = sum(len(m.group()) for m in re.finditer(r"((\w{1,3})\2{6,})", seq.upper())) / max(1, len(seq))
    return gc, low, simple
seen, rows = set(), []
loci = {}
for l in open(prof):
    f = l.rstrip("\n").split("\t")
    if f[0] == "gene" or f[9] != "deletion_larger_than_flanks" or not f[12] or int(f[12]) == 0:
        continue
    if f[4] in seen:
        continue
    seen.add(f[4])
    rows.append((f[0], f[4]))
# reference contig of each gene comes from verdicts via the profile coords; read verdicts
gene2c = {}
for l in open("out/verdicts.tsv"):
    f = l.split("\t")
    gene2c[f[0]] = f[2].split(":")[0]
print("gene\tbetulina_stretch\tlen\tGC\tlowercase\tsimple_repeat")
sizes = []
for g, coords in rows:
    c = gene2c[g]; s, e = map(int, coords.split("-"))
    gc, low, simple = stats(fetch(c, s, e))
    sizes.append((c, e - s))
    print(f"{g}\t{c}:{s}-{e}\t{e-s}\t{gc:.3f}\t{low:.3f}\t{simple:.4f}")
random.seed(1)
bg = []
for c, L in sizes * 4:
    st = random.randint(0, 5_000_000)
    bg.append(stats(fetch(c, st, st + L)))
print("background (random windows, same contigs/sizes): GC %.3f lowercase %.3f simple %.4f" % tuple(sum(x[i] for x in bg) / len(bg) for i in range(3)))
