#!/usr/bin/env python3
"""Re-label the weaker gap classes from read evidence (after recheck_reads.py).

Per h1 window (gap = N run x..y; tandem block aL..x / y..bR = merged TRF arrays >= 200 bp (gaps <= 50 bp) touching the gap, else
aL = x, bR = y):
  unique-anchored span: a read whose alignments (one strand, read order consistent) start <= aL - 500 and end
       >= bR + 500, with an alignment on each side of the N run (or one alignment straight through). Content = read
       bases between aL and bR minus the assembled array bases (= real gap content on that molecule; negative = the
       two flanks overlap, i.e. sequence written on both sides of the N run).
  hap: from recheck_reads.py events/summary (h1 = deposited hap1 base at >= 80 % of het sites, h2 = other allele)
  ramp: depth (h1-allele reads, h2-allele reads, all) at 15-35 kb vs 0-1 kb from the gap (depth.tsv)
New label (first that applies):
  both_read_through     >= 2 h1 and >= 2 h2 spanning reads (or >= 50 % of far depth spans when no het sites)
  other_hap_through     >= 2 h2 spans, <= 1 h1 span (deposited-haplotype reads fade -> its allele is unread)
  dep_hap_through       >= 2 h1 spans, <= 1 h2 span (the other haplotype's reads fade)
  through_unphased      >= 2 spans, haplotype unknown
  array_longer_than_reads  no span; a tandem block touches the gap and reads entering it from unique sequence end
                        inside array-like sequence after >= 10 kb of array, or the assembled block is >= 10 kb
  both_copies_drop      no span; all-read depth at the gap <= 40 % of the depth 15-35 kb away
  unspanned_other       no span, none of the above
Usage: recheck_label.py WIN.fa H1.paf H1.trf.dat RES_DIR TARGETS.tsv > relabel.tsv"""
import csv
import statistics
import sys
from collections import defaultdict

win, paf, trf, res, targets = sys.argv[1:6]
EDGE = 300
seqlen, xy = {}, {}
name = None
for l in open(win):
    if l.startswith(">"):
        name = l[1:].strip(); seqlen[name] = 0
        f = name.split("|")
        if f[1] == "h1":
            xy[name] = tuple(map(int, f[3].split("=")[1].split("-")))
    else:
        seqlen[name] += len(l.strip())
arr = defaultdict(list)
cur = None
for l in open(trf):
    if l.startswith("@"):
        cur = l[1:].strip()
    elif cur and l[:1].isdigit():
        f = l.split()
        if int(f[1]) - int(f[0]) + 1 >= 200:
            arr[cur].append([int(f[0]) - 1, int(f[1])])


def block(n):
    x, y = xy[n]
    m = []
    for a, b in sorted(arr[n]):
        if m and a <= m[-1][1] + 50:
            m[-1][1] = max(m[-1][1], b)
        else:
            m.append([a, b])
    aL, bR = x, y
    for a, b in m:
        if b >= x - 2000 and a < x:
            aL = min(aL, a)
        if a <= y + 2000 and b > y:
            bR = max(bR, b)
    return aL, bR


recs = defaultdict(lambda: defaultdict(list))
for l in open(paf):
    f = l.split("\t")
    if f[5] in xy:
        recs[f[5]][f[0]].append((int(f[7]), int(f[8]), int(f[2]), int(f[3]), f[4], int(f[1])))
hap = {}
for r in csv.DictReader(open(f"{res}/events.tsv"), delimiter="\t"):
    hap[(r["gap"], r["read"])] = r["hap"]
dep = defaultdict(list)
for r in csv.DictReader(open(f"{res}/depth.tsv"), delimiter="\t"):
    dep[r["gap"]].append((abs(int(r["dist_to_gap"])), int(r["all"]), int(r["h1"]), int(r["h2"])))
tand = defaultdict(dict)
for r in csv.DictReader(open(f"{res}/tandem.tsv"), delimiter="\t"):
    tand[r["gap"]][r["side"]] = r
cause = {r["gap"]: r["cause"] for r in csv.DictReader(open(targets), delimiter="\t")}


def mean(v):
    return sum(v) / len(v) if v else 0.0


print("gap\told_label\tarray_L\tarray_R\tn_het_reads\tspan_h1\tspan_h2\tspan_unph\tcontent_median\tcontent_list\t"
      "far_all\tnear_all\tfar_h1\tnear_h1\tfar_h2\tnear_h2\tarray_enter\tarray_end_inside\tarray_bp_max\tnew_label")
for n, (x, y) in sorted(xy.items()):
    g = n.split("|")[0]
    aL, bR = block(n)
    sp = defaultdict(list)
    for q, rr in recs[n].items():
        for st in "+-":
            r = sorted([v for v in rr if v[4] == st])
            if not r:
                continue
            left = [v for v in r if v[0] <= aL - 500 and v[1] >= min(x - EDGE, aL + 200)] if aL < x else \
                   [v for v in r if v[0] <= x - 1000 and v[1] >= x - EDGE]
            right = [v for v in r if v[1] >= bR + 500 and v[0] <= max(y + EDGE, bR - 200)] if bR > y else \
                    [v for v in r if v[1] >= y + 1000 and v[0] <= y + EDGE]
            thru = [v for v in r if v[0] <= min(aL - 500, x - 1000) and v[1] >= max(bR + 500, y + 1000)]
            if thru:
                v = thru[0]
                content = (v[3] - v[2]) - (v[1] - v[0]) + (y - x)
            elif left and right:
                a, b = left[0], right[-1]
                both = [v for v in left if v in right]
                if both:                          # one alignment satisfies both sides: straight through
                    v = both[0]
                    content = (v[3] - v[2]) - (v[1] - v[0]) + (y - x)
                    h = hap.get((g, q), "?")
                    sp[h].append(content)
                    break
                if a is b:
                    continue
                # read position of aL and bR
                if st == "+":
                    pa, pb = a[2] + (aL - a[0]), b[3] - (b[1] - bR)
                    if pb <= a[2]:
                        continue
                else:
                    pa, pb = a[3] - (aL - a[0]), b[2] + (b[1] - bR)
                    if pb >= a[3]:
                        continue
                content = abs(pb - pa) - ((bR - aL) - (y - x))
            else:
                continue
            h = hap.get((g, q), "?")
            sp[h].append(content)
            break
    D = dep[g]
    far = [d for d in D if 15000 <= d[0] < 35000]
    near = [d for d in D if d[0] < 1000]
    fa, na = mean([d[1] for d in far]), mean([d[1] for d in near])
    f1, n1 = mean([d[2] for d in far]), mean([d[2] for d in near])
    f2, n2 = mean([d[3] for d in far]), mean([d[3] for d in near])
    nh1, nh2, nu = len(sp["h1"]), len(sp["h2"]), len(sp["?"])
    allc = sorted(sp["h1"] + sp["h2"] + sp["?"])
    t = tand[g]
    enter = sum(int(t[s]["n_enter"]) for s in t)
    inside = sum(int(t[s]["in_array_end"]) for s in t)
    abmax = max([int(t[s]["array_bp_in_reads(max,median)"].split(",")[0]) for s in t
                 if t[s]["array_bp_in_reads(max,median)"]] or [0])
    if nh1 >= 2 and nh2 >= 2:
        lab = "both_read_through"
    elif nh2 >= 2 and nh1 <= 1:
        lab = "other_hap_through"
    elif nh1 >= 2 and nh2 <= 1:
        lab = "dep_hap_through"
    elif len(allc) >= 2 and (nh1 + nh2 == 0 and len(allc) >= 0.5 * fa):
        lab = "both_read_through"
    elif len(allc) >= 2:
        lab = "through_unphased"
    elif (aL < x - 0 or bR > y) and (abmax >= 10000 and inside >= 0.5 * max(1, enter) or (x - aL) >= 10000 or (bR - y) >= 10000):
        lab = "array_longer_than_reads"
    elif fa > 0 and na <= 0.4 * fa:
        lab = "both_copies_drop"
    else:
        lab = "unspanned_other"
    med = statistics.median(allc) if allc else ""
    print(f"{g}\t{cause.get(g)}\t{x - aL}\t{bR - y}\t-\t{nh1}\t{nh2}\t{nu}\t{med}\t{','.join(map(str, allc[:25]))}\t"
          f"{fa:.1f}\t{na:.1f}\t{f1:.1f}\t{n1:.1f}\t{f2:.1f}\t{n2:.1f}\t{enter}\t{inside}\t{abmax}\t{lab}")
