#!/usr/bin/env python3
"""One row per contig end of the betulina primary with measurements and a cause label (trizona gap_classify.py rules,
adapted to contig ends and ONT read-mapping depth instead of HiFi k-mers).

Depth: samtools bedcov -Q 20 over 1 kb bins (ends_prep.py); obs/exp = depth / (REF x two-copy expectation at the GC
of the end's outer 10 kb) (gc_curve.py on het bins). near = kb 0-1, mid = kb 2-4, far = median kb 10-29.
Sequence (outer NEAR = 2000 bp): telomere hexamer share, TRF tandem share, longest >= 90 % GA/TC tract, G4 motifs.
Context: het SNPs/kb in the outer 10 kb bin (Clair3); contig hemizygous (median 10 kb depth/REF <= 0.65 on contigs
>= 200 kb = X in this male); reads overhanging the end and where they continue (ends_reads.py); bp by which hap1 /
hap2 contigs continue past the end (end flank mapped to hap contigs); trizona hap1 position of the end and the nearest
trizona gap (<= 20 kb) with its cause, or a trizona sequence end (<= 50 kb).
Cause (first that applies):
  telomere        >= 30 % of the outer 2 kb is (TTAGGG)n / (CCCTAA)n -> real chromosome end
  tandem          >= 30 % tandem repeat (TRF)
  ga_tract        >= 200 bp GA/TC tract
  chrX            contig hemizygous (one copy in a male)
  high_copy       mid or far obs/exp >= 1.8 (collapsed repeat)
  unmappable      mid obs/exp < 0.3 (MAPQ >= 20 depth vanishes: repeat the reads cannot place)
  one_copy_fade   mid <= 0.75 while far >= 0.8 (coverage of one copy fades toward the end)
  one_copy_region far < 0.7 (whole flank one-copy: hemizygous segment / large het SV)
  other
Usage: ends_classify.py ENDS.fa ENDS.bedcov CURVE HET.tsv DEPTH10k.bed ENDS.trf.dat ENDS.reads.tsv QA_vs_H1.paf
       QA_vs_H2.paf ENDS_vs_T1.paf T1.fai TGAPS.tsv TCAUSES.tsv > end_causes.tsv"""
import re
import statistics
import sys
from collections import Counter, defaultdict

from calib import load_curve, read_fa, gc, med

(efa, bcov, curve, hetf, dep10, trf, rds, ph1, ph2, pt1, t1fai, tgaps, tcause) = sys.argv[1:14]
NEAR = 2000
expected, REF = load_curve(curve)
seqs = {k: v for k, v in read_fa(efa).items() if k.startswith("ptg")}
dep = defaultdict(dict)
for l in open(bcov):
    c, s, e, n, v = l.split()
    ctg, side, k = n.split("|")
    dep[f"{ctg}|{side}"][int(k)] = float(v) / (int(e) - int(s))
het = {}
for l in open(hetf):
    c, b, n = l.split()
    het[(c, int(b))] = int(n)
cd, clen = defaultdict(list), {}
for l in open(dep10):
    c, s, e, d = l.split()[:4]
    cd[c].append(float(d)); clen[c] = int(e)
hemi = {c for c, v in cd.items() if clen[c] >= 200000 and statistics.median(v) / REF <= 0.65}
tr, cur = defaultdict(list), None
for l in open(trf):
    if l.startswith("@"):
        cur = l[1:].split()[0]
    elif cur and l[:1].isdigit():
        f = l.split()
        tr[cur].append((int(f[0]) - 1, int(f[1])))
reads = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in list(open(rds))[1:]}


def outer(f):
    """PAF fields of an end-flank alignment -> (name, matches, target, outer point on target, bp beyond on target)"""
    q, ql, qs, qe, st, t, tl, ts, te, nm = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), \
        int(f[7]), int(f[8]), int(f[9])
    if q.endswith("|R"):  # outer side = flank end
        p = te + (ql - qe) if st == "+" else ts - (ql - qe)
        beyond = tl - p if st == "+" else p
    else:  # outer side = flank start
        p = ts - qs if st == "+" else te + qs
        beyond = p if st == "+" else tl - p
    return q, nm, t, p, beyond


def best_hits(paf):
    best = {}
    for l in open(paf):
        f = l.split("\t")
        if not f[0].startswith("ptg") or int(f[11]) < 5 or int(f[10]) < 2000:
            continue
        q, nm, t, p, beyond = outer(f)
        if q not in best or nm > best[q][0]:
            best[q] = (nm, t, p, beyond)
    return best


h1, h2, t1 = best_hits(ph1), best_hits(ph2), best_hits(pt1)
tl = {l.split()[0]: int(l.split()[1]) for l in open(t1fai)}
tg = defaultdict(list)
for l in list(open(tgaps))[1:]:
    f = l.split("\t")
    tg[f[1]].append((int(f[2]), int(f[3]), f[0]))
tc = {l.split("\t")[0]: l.rstrip("\n").split("\t")[-1] for l in list(open(tcause))[1:]}


def tri_ctx(e):
    if e not in t1:
        return "", "", "unmapped"
    _, t, p, _ = t1[e]
    best = None
    for a, b, g in tg.get(t, []):
        d = 0 if a <= p <= b else min(abs(p - a), abs(p - b))
        if d <= 20000 and (best is None or d < best[0]):
            best = (d, g)
    if best:
        return t, p, f"gap:{best[1]}:{tc.get(best[1], '?')}"
    if min(p, tl[t] - p) <= 50000:
        return t, p, "seq_end"
    return t, p, "contiguous"


def feats(e):
    s = seqs[e]
    R = e.endswith("|R")
    a = max(0, len(s) - NEAR) if R else 0
    x = s[a:a + NEAR]
    telo = (len(re.findall("TTAGGG", x)) + len(re.findall("CCCTAA", x))) * 6 / max(1, len(x))
    cov = [0] * len(x)
    for s0, e0 in tr[e]:
        for i in range(max(a, s0), min(a + len(x), e0)):
            cov[i - a] = 1
    best = 0
    for alph in ("AG", "CT"):
        i = 0
        while i < len(x):
            j, bad = i, 0
            while j < len(x) and (bad + (x[j] not in alph)) <= 0.1 * (j - i + 1):
                bad += x[j] not in alph; j += 1
            best = max(best, j - i); i = max(j, i + 1)
    g4 = len(re.findall(r"(?=G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})", x)) + \
        len(re.findall(r"(?=C{3,}\w{1,7}C{3,}\w{1,7}C{3,}\w{1,7}C{3,})", x))
    g10 = gc(s[-10000:] if R else s[:10000])
    return telo, sum(cov) / max(1, len(cov)), best, g4, g10


cols = ["end", "ctg_len", "gc10k", "near_oe", "mid_oe", "far_oe", "het_per_kb", "hemi", "telo", "tandem", "ga_run", "g4",
        "n_reach", "n_overhang", "read_partner", "n_partner", "hap1_beyond", "hap2_beyond", "tri_chrom", "tri_pos",
        "tri_context", "cause"]
print("\t".join(cols))
cnt = Counter()
f2 = lambda v: "" if v != v else f"{v:.2f}"
for e in seqs:
    c, side = e.split("|")
    r = reads.get(e, [e, "0", "", "", "-", "0", "0"])
    L = int(r[1])
    telo, tand, ga, g4, g10 = feats(e)
    ex = expected(round(g10)) * REF
    d = dep[e]
    oe = lambda ks: [d[k] / ex for k in ks if k in d]
    near, mid, far = med(oe([0, 1])), med(oe([2, 3, 4])), med(oe(range(10, 30)))
    hb = 0 if side == "L" else ((L - 1) // 10000) * 10000
    hk = het.get((c, hb), 0) / 10
    b1 = h1[e][3] if e in h1 else ""
    b2 = h2[e][3] if e in h2 else ""
    tchr, tpos, tctx = tri_ctx(e)
    if telo >= 0.3:
        cause = "telomere"
    elif tand >= 0.3:
        cause = "tandem"
    elif ga >= 200:
        cause = "ga_tract"
    elif c in hemi:
        cause = "chrX"
    elif (mid == mid and mid >= 1.8) or (far == far and far >= 1.8):
        cause = "high_copy"
    elif mid == mid and mid < 0.3:
        cause = "unmappable"
    elif mid == mid and mid <= 0.75 and far == far and far >= 0.8:
        cause = "one_copy_fade"
    elif far == far and far < 0.7:
        cause = "one_copy_region"
    else:
        cause = "other"
    cnt[cause] += 1
    row = [e, L, f"{g10:.1f}", f2(near), f2(mid), f2(far), f"{hk:.2f}", int(c in hemi), f"{telo:.2f}", f"{tand:.2f}", ga,
           g4, r[2], r[3], r[4], r[5], b1, b2, tchr, tpos, tctx, cause]
    print("\t".join(map(str, row)))
for k, v in cnt.most_common():
    print(f"# {k}\t{v}", file=sys.stderr)
