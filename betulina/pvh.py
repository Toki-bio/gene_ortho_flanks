#!/usr/bin/env python3
"""Primary vs hap1 / hap2: which primary 10 kb bins are carried by hap1, hap2, both or neither, as a function of
heterozygosity and GC-calibrated read depth (trizona mechanism A: a homozygous stretch whose coverage falls below
hifiasm's homozygous threshold is treated as haplotype-specific and put into ONE haplotype).
Query may be chunk_fa.py chunks (name seq:offset). A bin is 'in hapN' when >= 50 % of it is covered by query-side intervals of primary alignments (tp:A:P, >= 5 kb) of the
primary contigs to hapN (minimap2 -x asm5 HAP PRIMARY).
Autosomal bins only (contig >= 1 Mb, not hemizygous). het = >= 1 het SNP/kb, hetfree = < 0.3/kb.
Usage: pvh.py P_vs_H1.paf P_vs_H2.paf DEPTH10k.bed GC.bed HET.tsv CURVE"""
import statistics
import sys
from collections import defaultdict

from calib import load_curve

p1, p2, dep10, gcf, hetf, curve = sys.argv[1:7]
BIN = 10000
expected, REF = load_curve(curve)


def covered(paf):
    cov = defaultdict(int)
    for l in open(paf):
        f = l.split("\t")
        if "tp:A:P" not in l or int(f[3]) - int(f[2]) < 5000:
            continue
        q, a, e = f[0], int(f[2]), int(f[3])
        if ":" in q:  # chunk_fa.py chunk name seq:offset
            q, off = q.rsplit(":", 1); a += int(off); e += int(off)
        while a < e:
            b = a // BIN * BIN; k = min(e, b + BIN) - a
            cov[(q, b)] += k; a += k
    return cov


c1, c2 = covered(p1), covered(p2)
gc = {}
for l in open(gcf):
    c, s, e, g, n = l.split()
    if float(g) >= 0 and int(n) >= 0.8 * (int(e) - int(s)):
        gc[(c, int(s))] = (round(float(g)), int(e) - int(s))
het = {(l.split()[0], int(l.split()[1])): int(l.split()[2]) for l in open(hetf)}
rows, cd, clen = [], defaultdict(list), {}
for l in open(dep10):
    c, s, e, d = l.split()[:4]
    cd[c].append(float(d)); clen[c] = int(e)
    rows.append((c, int(s), float(d)))
hemi = {c for c, v in cd.items() if clen[c] >= 200000 and statistics.median(v) / REF <= 0.65}
tab = defaultdict(lambda: defaultdict(int))
gtab = defaultdict(lambda: defaultdict(int))
runs, cur = [], None
for c, s, d in rows:
    k = (c, s)
    if c in hemi or clen[c] < 1000000 or k not in gc:
        continue
    g, w = gc[k]
    oe = d / (REF * expected(g))
    h = het.get(k, 0) / w * 1000
    st = ("both" if c1[k] >= w / 2 and c2[k] >= w / 2 else "h1_only" if c1[k] >= w / 2 else
          "h2_only" if c2[k] >= w / 2 else "neither")
    hz = "het" if h >= 1 else "hetfree" if h < 0.3 else "mid"
    ob = min(int(oe / 0.1) / 10, 1.6)
    tab[(hz, ob)][st] += 1
    if hz == "hetfree" and 0.3 <= oe <= 1.6:
        gtab[min(max(g // 2 * 2, 34), 56)][st] += 1
    one = st in ("h1_only", "h2_only")
    if one and cur and cur[0] == c and cur[2] == s and cur[3] == st:
        cur[2] = s + w; cur[4].append(oe); cur[5].append(h)
    else:
        if cur:
            runs.append(cur)
        cur = [c, s, s + w, st, [oe], [h]] if one else None
if cur:
    runs.append(cur)
STS = ["both", "h1_only", "h2_only", "neither"]
print("# bins by het state and depth obs/exp (two-copy = 1): share carried by both haps / one hap / neither")
print("het\toe_bin\tn\t" + "\t".join(STS) + "\tone_hap_share")
for (hz, ob) in sorted(tab):
    t = tab[(hz, ob)]; n = sum(t.values())
    if n >= 50:
        print(f"{hz}\t{ob:.1f}\t{n}\t" + "\t".join(str(t[s]) for s in STS) + f"\t{(t['h1_only'] + t['h2_only']) / n:.3f}")
print("# het-free bins (oe 0.3-1.6) by GC: share in one hap only")
for g in sorted(gtab):
    t = gtab[g]; n = sum(t.values())
    print(f"GC{g}\t{n}\t" + "\t".join(str(t[s]) for s in STS) + f"\t{(t['h1_only'] + t['h2_only']) / n:.3f}")
big = [r for r in runs if r[2] - r[1] >= 30000]
print(f"# one-hap-only runs >= 30 kb: {len(big)}, total {sum(r[2] - r[1] for r in big) / 1e6:.1f} Mb; "
      f"h1_only {sum(1 for r in big if r[3] == 'h1_only')}, h2_only {sum(1 for r in big if r[3] == 'h2_only')}")
hf = [r for r in big if statistics.median(r[5]) < 0.3]
print(f"#   of these het-free (median < 0.3 het/kb): {len(hf)}, {sum(r[2] - r[1] for r in hf) / 1e6:.1f} Mb, "
      f"median depth obs/exp {statistics.median([statistics.median(r[4]) for r in hf]) if hf else float('nan'):.2f}; "
      f"het runs median obs/exp "
      f"{statistics.median([statistics.median(r[4]) for r in big if statistics.median(r[5]) >= 1]) if any(statistics.median(r[5]) >= 1 for r in big) else float('nan'):.2f}")
