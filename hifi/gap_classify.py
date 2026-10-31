#!/usr/bin/env python3
"""One row per gap of the gapped assembly with everything measured, and a cause label.

Inputs: gapfill.py gaps.tsv (class), gap_flank_cn.py table (fill call, fill and flank copy number), the flank sequences
(gap_flank_cn.py prep output), TRF -ngs on those flanks, and the X chromosome name.
Per flank, the NEAR bp next to the gap: GC, tandem share (TRF), longest purine (A/G) or pyrimidine (C/T) run with
>= 90 % purity (GA/TC tracts are a documented HiFi dropout), and the flank copy number.
Cause labels (first that applies):
  hom_one_hap   a fill exists in another assembly, its copy number is two-copy  -> homozygous stretch put into one
                haplotype (GC-depressed coverage below hifiasm's homozygous threshold)
  one_copy_fill a fill exists, copy number one-copy                         -> only one copy carries it
  join          the flanks meet directly in another assembly (the Ns hide nothing)
  tandem        unspanned, a flank's NEAR bp is >= 30 % tandem repeat
  ga_tract      unspanned, a >= 200 bp GA/TC tract in a flank's NEAR bp
  chrX          unspanned, on the X (one copy in a male: half coverage)
  high_copy     unspanned, flank copy number >= 1.8 (collapsed repeat)
  other         none of the above
Usage: gap_classify.py GAPS.tsv FLANK_CN.tsv FLANKS.fa FLANKS.trf.dat X_NAME > gap_causes.tsv"""
import re
import sys
from collections import Counter, defaultdict

gtsv, fcn, ffa, trf, xname = sys.argv[1:6]
NEAR = 2000
seqs, name = {}, None
for l in open(ffa):
    if l.startswith(">"):
        name = l[1:].strip(); seqs[name] = ""
    else:
        seqs[name] += l.strip().upper()
tr = defaultdict(list); cur = None
for l in open(trf):
    if l.startswith("@"):
        cur = l[1:].split()[0]
    elif cur and l[:1].isdigit():
        f = l.split()
        tr[cur].append((int(f[0]) - 1, int(f[1])))


def near_feats(n):
    s = seqs.get(n)
    if not s:
        return None
    L = len(s)
    a, b = (L - NEAR, L) if n.endswith("|L") else (0, NEAR)
    x = s[a:b]
    gc = 100 * (x.count("G") + x.count("C")) / max(1, sum(x.count(y) for y in "ACGT"))
    cov = [0] * (b - a)
    for s0, e0 in tr[n]:
        for i in range(max(a, s0), min(b, e0)):
            cov[i - a] = 1
    best = 0
    for alph in ("AG", "CT"):
        # longest window where >= 90 % of bases are in alph (greedy scan with a 50 bp seed)
        run = 0
        for m in re.finditer(f"[{alph}]{{50,}}", x):
            run = max(run, len(m.group()))
        ext = 0
        i = 0
        while i < len(x):
            j, bad = i, 0
            while j < len(x) and (bad + (x[j] not in alph)) <= 0.1 * (j - i + 1):
                bad += x[j] not in alph; j += 1
            ext = max(ext, j - i); i = max(j, i + 1)
        best = max(best, run, ext)
    return gc, sum(cov) / max(1, len(cov)), best


cn = {}
for l in list(open(fcn))[1:]:
    f = l.rstrip("\n").split("\t")
    cn[f[0]] = f
chrom = {}
for l in list(open(gtsv))[1:]:
    f = l.split("\t")
    chrom[f[0]] = f[1]
print("gap\tchrom\tclass\tfill_call\tfill_oe\tflankL_oe\tflankR_oe\tnearGC_L\tnearGC_R\ttandem_L\ttandem_R\tgaRun_L\tgaRun_R\tcause")
counts = Counter()
for g, c in chrom.items():
    f = cn.get(g, [g, "", "-", "", "", ""])
    cls, fcall = f[1], f[2]
    fl = [float(x) for x in f[4:6] if x]
    ft = {s: near_feats(f"{g}|{s}") for s in "LR"}
    tand = max((v[1] for v in ft.values() if v), default=0)
    ga = max((v[2] for v in ft.values() if v), default=0)
    if cls == "filled" and fcall == "two_copies":
        cause = "hom_one_hap"
    elif cls == "filled" and fcall == "one_copy":
        cause = "one_copy_fill"
    elif cls == "filled":
        cause = f"fill_{fcall.replace('-', 'nocall')}"
    elif cls == "joined":
        cause = "join"
    elif tand >= 0.3:
        cause = "tandem"
    elif ga >= 200:
        cause = "ga_tract"
    elif c == xname:
        cause = "chrX"
    elif fl and max(fl) >= 1.8:
        cause = "high_copy"
    else:
        cause = "other"
    counts[(cls, cause)] += 1
    val = lambda s, i: f"{ft[s][i]:.2f}" if ft[s] else ""
    print(f"{g}\t{c}\t{cls}\t{fcall}\t{f[3]}\t{f[4]}\t{f[5]}\t{val('L', 0)}\t{val('R', 0)}\t{val('L', 1)}\t{val('R', 1)}\t"
          f"{ft['L'][2] if ft['L'] else ''}\t{ft['R'][2] if ft['R'] else ''}\t{cause}")
for (cls, cause), n in sorted(counts.items()):
    print(f"# {cls}\t{cause}\t{n}", file=sys.stderr)
