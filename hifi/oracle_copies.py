#!/usr/bin/env python3
"""How many copies of each traced region does an independent HiCanu assembly hold?

A region is a stretch hifiasm wrote in one haplotype only. One HiCanu contig matching it proves nothing, because that
contig may just be the copy hifiasm kept. So each window is cut into PIECE-bp pieces and every piece is mapped alone
(a whole-contig alignment hides the other parent's contig: secondary chains under 0.8 of a long primary are dropped).
Per piece count the contigs carrying it (aligned >= half the piece, identity >= MIN_ID) and add their read coverage
(tigInfo): HiCanu merges identical parental copies into one contig, so one contig at ~2 x HAP_COV is two copies.
Copy estimate of a piece = sum of contig coverages / HAP_COV, rounded. Contig coverage is a whole-contig mean, so short
or repeat-rich contigs make it rough. The windows labelled CTRL are ordinary two-copy regions: they must read 2.

  oracle_copies.py pieces WINDOWS.bed HAP1.fa HAP2.fa > pieces.fa
  minimap2 -c -x asm20 -p 0.5 -N 20 CONTIGS.fa pieces.fa > pieces.paf
  oracle_copies.py count WINDOWS.bed TIGINFO HAP_COV pieces.paf > copies.tsv"""
import sys
from collections import defaultdict

PIECE = 20000
MIN_ID = 0.98
mode, win = sys.argv[1:3]
windows = [l.split() for l in open(win)]

if mode == "pieces":
    fa = {"hap1": sys.argv[3], "hap2": sys.argv[4]}
    for hap in fa:
        need = defaultdict(list)
        for h, sq, s, e, label in windows:
            if h == hap:
                need[sq].append((int(s), int(e), label))
        name, buf = None, []

        def flush():
            if name in need:
                seq = "".join(buf)
                for s, e, label in need[name]:
                    for b in range(s, e, PIECE):
                        print(f">{label}|{hap}|{b}\n{seq[b:min(b + PIECE, e)]}")

        for l in open(fa[hap]):
            if l[0] == ">":
                flush()
                name, buf = l[1:].split()[0], []
            elif name in need:
                buf.append(l.strip())
        flush()
    sys.exit()

tiginfo, hap_cov, paf = sys.argv[3], float(sys.argv[4]), sys.argv[5]
cov = {}
for l in open(tiginfo):
    if l[0] != "#":
        f = l.split("\t")
        cov[int(f[0])] = float(f[2])
hits = defaultdict(lambda: defaultdict(int))     # piece -> contig -> aligned bp
plen = {}
for l in open(paf):
    f = l.split("\t")
    plen[f[0]] = int(f[1])
    if int(f[9]) / int(f[10]) >= MIN_ID:
        hits[f[0]][f[5]] += int(f[3]) - int(f[2])
print("label\thap\tseq\tstart\tend\tpieces\tshare0\tshare1\tshare2\tshare3+\tcontigs\tprofile")
for hap, sq, s, e, label in windows:
    s, e = int(s), int(e)
    est = []
    used = set()
    for b in range(s, e, PIECE):
        p = f"{label}|{hap}|{b}"
        keep = [t for t, v in hits[p].items() if v >= min(b + PIECE, e) - b >> 1]
        used.update(keep)
        est.append(round(sum(cov.get(int(t[3:]), 0) for t in keep) / hap_cov))
    n = len(est)
    sh = [sum(1 for x in est if x == k) / n for k in (0, 1, 2)]
    sh.append(sum(1 for x in est if x >= 3) / n)
    prof = "".join(str(min(x, 9)) for x in est)
    print(f"{label}\t{hap}\t{sq}\t{s}\t{e}\t{n}\t" + "\t".join(f"{x:.2f}" for x in sh) + f"\t{len(used)}\t{prof}")
