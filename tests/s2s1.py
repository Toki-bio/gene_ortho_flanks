#!/usr/bin/env python3
"""Distribution of s2/s1 (minimap2: best secondary / primary chaining score) of the
best primary piece per flank, by gene group. Usage: s2s1.py OUT_DIR TEST_GENES FI_GENES"""
import collections
import re
import sys

out, test_f, fi_f = sys.argv[1:4]
test = set(l.split()[0] for l in open(test_f) if l.strip())
fi = set(l.split()[0] for l in open(fi_f) if l.strip())
ig = re.compile(r"^(IG[HKL]V|IG[HKL][DJ]|TR[ABGD]V|TR[ABGD]J)")
best = {}
for l in open(out + "/flanks.paf"):
    if "tp:A:P" not in l:
        continue
    f = l.rstrip("\n").split("\t")
    tags = dict(t.split(":", 2)[::2] for t in f[12:])
    s1, s2 = int(tags.get("s1", 0)), int(tags.get("s2", 0))
    if f[0] not in best or s1 > best[f[0]][0]:
        best[f[0]] = (s1, s2)
bins = (0, 0.5, 0.8, 0.9, 0.95, 1.0, float("inf"))  # s2 can exceed s1 after DP rescoring
tab = collections.defaultdict(collections.Counter)
for n, (s1, s2) in best.items():
    g, side, fl = n.split("|")
    grp = "FI" if g in fi else ("IG/TR" if ig.match(g) else "other_LM") if g in test else None
    if grp is None:
        continue
    r = s2 / s1 if s1 else 0
    b = next(f"{bins[i]}-{bins[i+1]}" for i in range(len(bins) - 1) if r < bins[i + 1])
    tab[(grp, fl)][b] += 1
labels = [f"{bins[i]}-{bins[i+1]}" for i in range(len(bins) - 1)]
print("group\tflank\t" + "\t".join(labels))
for k in sorted(tab):
    print(*k, *(tab[k][b] for b in labels), sep="\t")
