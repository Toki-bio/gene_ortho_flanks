#!/usr/bin/env python3
"""What are the 'fill is a repeat' fills (fill_cn.py call `repeat`: median fill k-mer CN >= 2x the two-copy level)?

Two steps.
  prep:  recheck_fillrep.py prep ALL_FILLS.fa FILL_CN.tsv CAUSES.tsv CAUSE OUT
         -> OUT/core.fa (one fill per gap, the 1 kb source-contig ends removed; dep2 preferred, else the first source
            whose call is `repeat`), OUT/core.ids
  table: recheck_fillrep.py table OUT ALL_FILLS.wig CALIB TRF.dat DFAM.m8 DFAM_META SELF.m8 GENOME.paf [PEAK]
         DFAM.m8 / SELF.m8: blastn -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore'
         GENOME.paf: minimap2 -c of core.fa to the deposited hap1 (the fill is absent from hap1, so every hit is
         another copy)
Per fill (position-wise k-mer CN r = count / PEAK / two-copy expectation at the fill's GC):
  hi     share of positions with r >= 2 (repeat k-mers), uniq share with 0.3 <= r < 1.5
  r_uniq median r over the uniq positions (~1 two copies, ~0.5 one copy)
  trf    share covered by TRF arrays (period, copies, % match of the longest array)
  dfam   share covered by Dfam curated consensus hits (e <= 1e-5), the top family, its % identity and the share of
         the consensus it covers
  self   share of the fill that matches elsewhere in the fill (internal duplication), excluding the diagonal
  gcov   share of the fill covered by hap1 alignments >= 1 kb at >= 90 % identity, number of such loci, best identity
Label (first that applies): tandem (trf >= 0.5), TE:<family> (dfam >= 0.5), segdup (gcov >= 0.5),
  mixed_repeat (hi >= 0.5), mostly_unique (hi < 0.5: the repeat call came from a minority of the fill)"""
import statistics
import sys
from collections import defaultdict

FO = 1000


def read_fa(path):
    seqs, name, buf = {}, None, []
    for l in list(open(path)) + [">"]:
        if l.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name, buf = (l[1:].split()[0] if len(l.strip()) > 1 else None), []
        else:
            buf.append(l.strip())
    return seqs


def cover(ivs, L):
    m = [0] * L
    for a, b in ivs:
        for i in range(max(0, a), min(L, b)):
            m[i] = 1
    return sum(m) / max(1, L)


def prep(argv):
    fa, fcn, causes, cause, out = argv
    want = {l.split("\t")[0] for l in open(causes) if l.rstrip("\n").split("\t")[-1] == cause}
    call = {}
    for l in list(open(fcn))[1:]:
        f = l.rstrip("\n").split("\t")
        call[(f[0], f[1])] = f[-1]
    seqs = read_fa(fa)
    pick = {}
    for n in seqs:
        g, src = n.split("|")[:2]
        if g not in want:
            continue
        if g not in pick or (src == "dep2" and pick[g].split("|")[1] != "dep2"):
            if call.get((g, src)) == "repeat" or g not in pick:
                pick[g] = n
    with open(f"{out}/core.fa", "w") as fh, open(f"{out}/core.ids", "w") as ih:
        for g in sorted(pick):
            s = seqs[pick[g]][FO:-FO]
            fh.write(f">{g}\n{s}\n"); ih.write(f"{g}\t{pick[g]}\n")
    print(f"{len(pick)} fills of {len(want)} gaps", file=sys.stderr)


def table(argv):
    out, wig, calib, trf, dfam, meta, selfm8, gpaf = argv[:8]
    PEAK = float(argv[8]) if len(argv) > 8 else 26
    exp = {}
    for x in open(calib).readline().split(":", 1)[1].split():
        g, v = x.split(":"); exp[int(g)] = float(v)
    core = read_fa(f"{out}/core.fa")
    src = dict(l.rstrip("\n").split("\t") for l in open(f"{out}/core.ids"))
    back = {v: k for k, v in src.items()}
    cnt = defaultdict(dict)
    chrom = None
    for l in open(wig):
        if l.startswith("variableStep"):
            chrom = back.get(l.split("chrom=")[1].split()[0]); continue
        if chrom:
            p, c = l.split()
            p = int(p) - 1 - FO
            if 0 <= p < len(core[chrom]):
                cnt[chrom][p] = int(c)
    tr = defaultdict(list)
    cur = None
    for l in open(trf):
        if l.startswith("@"):
            cur = l[1:].split()[0]
        elif cur and l[:1].isdigit():
            f = l.split()
            tr[cur].append((int(f[0]) - 1, int(f[1]), int(f[2]), float(f[3]), int(f[5])))
    name = {}
    for l in open(meta):
        f = l.split("\t")
        name[f[0]] = f[1]
    dh = defaultdict(list)
    for l in open(dfam):
        f = l.split("\t")
        if float(f[8]) <= 1e-5:
            dh[f[0]].append((int(f[4]) - 1, int(f[5]), f[1], float(f[2]), min(int(f[6]), int(f[7])), max(int(f[6]), int(f[7]))))
    sh = defaultdict(list)
    for l in open(selfm8):
        f = l.split("\t")
        if f[0] != f[1]:
            continue
        qs, qe, ss, se = int(f[4]), int(f[5]), int(f[6]), int(f[7])
        if qs == ss and qe == se:
            continue
        if int(f[3]) >= 100:
            sh[f[0]].append((qs - 1, qe))
    gh = defaultdict(list)
    for l in open(gpaf):
        f = l.split("\t")
        nm, bl = int(f[9]), int(f[10])
        if bl >= 1000 and nm / bl >= 0.90:
            gh[f[0]].append((int(f[2]), int(f[3]), f[5], int(f[7]), nm / bl))
    dlen = {}
    print("gap\tsource\tfill_len\tGC\thi\tuniq\tr_uniq\ttrf\ttrf_top(period,copies,%match,len)\tdfam\tdfam_top(family,%id,cons_cov)"
          "\tself\tgcov\tgloci\tgbest_id\tlabel")
    for g, s in core.items():
        L = len(s)
        acgt = sum(s.count(x) for x in "ACGT")
        gc = 100 * (s.count("G") + s.count("C")) / max(1, acgt)
        e = exp[min(exp, key=lambda k: abs(k - gc))]
        rs = [c / PEAK / e for c in cnt[g].values() if c > 0]
        hi = sum(r >= 2 for r in rs) / max(1, len(rs))
        un = [r for r in rs if 0.3 <= r < 1.5]
        ru = statistics.median(un) if un else float("nan")
        t = tr[g]
        tcov = cover([(a, b) for a, b, *_ in t], L)
        top = max(t, key=lambda x: x[1] - x[0], default=None)
        ttxt = f"{top[2]},{top[3]:.0f},{top[4]},{top[1] - top[0]}" if top else ""
        d = dh[g]
        dcov = cover([(a, b) for a, b, *_ in d], L)
        fam = defaultdict(int)
        for a, b, sid, pid, cs, ce in d:
            fam[sid] += b - a
        dtxt = ""
        if fam:
            sid = max(fam, key=fam.get)
            hs = [x for x in d if x[2] == sid]
            pid = sum((x[1] - x[0]) * x[3] for x in hs) / sum(x[1] - x[0] for x in hs)
            dtxt = f"{name.get(sid, sid)},{pid:.1f},{fam[sid]}bp"
        scov = cover(sh[g], L)
        gcv = cover([(a, b) for a, b, *_ in gh[g]], L)
        loci = len({(c, p // 100000) for _, _, c, p, _ in gh[g]})
        gb = max((x[4] for x in gh[g]), default=0)
        if tcov >= 0.5:
            lab = "tandem"
        elif dcov >= 0.5:
            lab = "TE:" + (dtxt.split(",")[0] if dtxt else "?")
        elif gcv >= 0.5:
            lab = "segdup"
        elif hi >= 0.5:
            lab = "mixed_repeat"
        else:
            lab = "mostly_unique"
        print(f"{g}\t{src[g].split('|')[1]}\t{L}\t{gc:.1f}\t{hi:.2f}\t{len(un) / max(1, len(rs)):.2f}\t{ru:.2f}\t{tcov:.2f}\t{ttxt}"
              f"\t{dcov:.2f}\t{dtxt}\t{scov:.2f}\t{gcv:.2f}\t{loci}\t{gb:.3f}\t{lab}")


if __name__ == "__main__":
    {"prep": prep, "table": table}[sys.argv[1]](sys.argv[2:])
