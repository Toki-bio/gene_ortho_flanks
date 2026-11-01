#!/usr/bin/env python3
"""Trizona problem sites (and random trizona sites) on betulina: does the betulina assembly break there, and is the
betulina ONT depth reduced there?

prep:    site_depth.py prep SITES_P.tsv PRIMARY.fa OUT.bed
         For sites 'spanned' in the primary: centre bins (1 kb) over [min point - 1 kb, max point + 1 kb] named
         site|C|i, and 1 kb flank bins 1..9 kb outward on both sides named site|F|k; GC of each bin in column 5.
summary: site_depth.py summary OUT.bed BEDCOV CURVE CAUSES.tsv SITES_P.tsv SITES_H1.tsv SITES_H2.tsv
         BEDCOV = samtools bedcov -Q 20 OUT.bed BAM (6th column = summed depth).
         Per cause class (trizona gap_causes.tsv; rnd* sites = 'random'): sites, both flanks placed, share spanned in
         primary / hap1 / hap2, share split at a contig end, centre min depth obs/exp (median, share <= 0.7,
         share <= 0.3) and the flank ramp (median obs/exp at 1, 3, 5, 9 kb)."""
import statistics
import sys
from collections import Counter, defaultdict

from calib import load_curve, read_fa, gc, med

if sys.argv[1] == "prep":
    sites, fa, out = sys.argv[2:5]
    rows = [l.rstrip("\n").split("\t") for l in list(open(sites))[1:]]
    rows = [r for r in rows if r[1] == "spanned"]
    need = {r[3] for r in rows}
    seqs = {k: v for k, v in read_fa(fa).items() if k in need}
    with open(out, "w") as fh:
        for r in rows:
            c, a, b = r[3], int(r[4]), int(r[7])
            lo, hi = max(0, min(a, b) - 1000), min(len(seqs[c]), max(a, b) + 1000)
            bins = [(f"{r[0]}|C|{i}", x, min(hi, x + 1000)) for i, x in enumerate(range(lo, hi, 1000))]
            for k in range(1, 10):
                bins.append((f"{r[0]}|F|{k}", lo - k * 1000, lo - (k - 1) * 1000))
                bins.append((f"{r[0]}|F|{k}", hi + (k - 1) * 1000, hi + k * 1000))
            for n, s, e in bins:
                if s >= 0 and e <= len(seqs[c]) and e - s >= 200:
                    fh.write(f"{c}\t{s}\t{e}\t{n}\t{gc(seqs[c][s:e]):.1f}\n")
    sys.exit(0)

bed, bcov, curve, causes, sp, s1, s2 = sys.argv[2:9]
expected, REF = load_curve(curve)
cause = {l.split("\t")[0]: l.rstrip("\n").split("\t")[-1] for l in list(open(causes))[1:]}
grp = lambda s: "random" if s.startswith("rnd") else cause.get(s, "?")
st, dist = {}, {}
for lab, f in (("P", sp), ("H1", s1), ("H2", s2)):
    for l in list(open(f))[1:]:
        r = l.split("\t")
        st.setdefault(r[0], {})[lab] = r[1]
        if lab == "P" and r[1] == "spanned":
            dist[r[0]] = int(r[2])
cen, ramp = defaultdict(list), defaultdict(lambda: defaultdict(list))
for l in open(bcov):
    c, s, e, n, g, v = l.split()
    oe = float(v) / (int(e) - int(s)) / (REF * expected(round(float(g))))
    site, kind, k = n.split("|")
    if kind == "C":
        cen[site].append(oe)
    else:
        ramp[site][int(k)].append(oe)
G = defaultdict(list)
for s in st:
    G[grp(s)].append(s)
print("class\tn\tplaced_P\tspan_P\tspan_H1\tspan_H2\tspan_any\tsplitEnd_P\tcen_min_oe_med\tcen_le0.7\tcen_le0.3\t"
      "cen_med_oe\tdist_med\tramp1\tramp3\tramp5\tramp9")
order = sorted(G, key=lambda g: (g == "random", -len(G[g])))
for g in order + ["ALL_GAPS"]:
    ss = [s for s in st if not s.startswith("rnd")] if g == "ALL_GAPS" else G[g]
    n = len(ss)
    placed = [s for s in ss if st[s].get("P") in ("spanned", "split_end", "split_other")]
    sh = lambda lab, val: sum(1 for s in placed if st[s].get(lab) == val) / max(1, len(placed))
    anysp = sum(1 for s in placed if "spanned" in st[s].values()) / max(1, len(placed))
    mins = [min(cen[s]) for s in placed if cen.get(s)]
    rk = lambda k: med([statistics.median(ramp[s][k]) for s in placed if ramp[s].get(k)])
    f = lambda v: f"{v:.2f}"
    print(f"{g}\t{n}\t{len(placed)}\t{f(sh('P', 'spanned'))}\t{f(sh('H1', 'spanned'))}\t{f(sh('H2', 'spanned'))}\t"
          f"{f(anysp)}\t{f(sh('P', 'split_end'))}\t{f(med(mins))}\t"
          f"{f(sum(v <= 0.7 for v in mins) / max(1, len(mins)))}\t{f(sum(v <= 0.3 for v in mins) / max(1, len(mins)))}\t"
          f"{f(med([statistics.median(cen[s]) for s in placed if cen.get(s)]))}\t"
          f"{med([dist[s] for s in placed if s in dist])}\t"
          f"{f(rk(1))}\t{f(rk(3))}\t{f(rk(5))}\t{f(rk(9))}")
