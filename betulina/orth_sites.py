#!/usr/bin/env python3
"""Where do the two flanks of a site (gap or random point) of one genome land in another assembly?
Input: PAF of the flank FASTA (names SITE|L = left flank ending at the site, SITE|R = right flank starting at it;
gap_flank_cn.py prep) against the assembly (minimap2 -x asm20), and the assembly .fai.
Per flank the best alignment (most matches) with MAPQ >= MINQ, >= MINALN bp, reaching within SLACK bp of the site side
of the flank. Its site point on the target = alignment edge extended by the unaligned flank bases on the site side.
Status per site:
  spanned     both flanks on one contig, same strand, in order; dist = bp between the two site points (target bases
              standing in for the site; ~0 for a simple point; can be large if the target has an insertion)
  split_end   both found, different contigs (or out of order), and >= 1 site point within END bp of its contig end
              -> the target assembly breaks here
  split_other both found, different contigs / order, neither at a contig end (rearrangement or paralog)
  left_only / right_only / none
Output: site  status  dist  contigL  pointL  strandL  contigR  pointR  strandR  endL  endR
Usage: orth_sites.py FLANKS.paf TARGET.fa.fai SITES_LIST [MAXDIST=200000] > sites.tsv
SITES_LIST: file whose first column is the site name (header line skipped)."""
import sys

paf, fai, sl = sys.argv[1:4]
MAXDIST = int(sys.argv[4]) if len(sys.argv) > 4 else 200000
MINQ, MINALN, SLACK, END = 5, 2000, 3000, 20000
ln = {l.split()[0]: int(l.split()[1]) for l in open(fai)}
best = {}
for l in open(paf):
    f = l.split("\t")
    q, ql, qs, qe, st, t, tl, ts, te, nm, al, mq = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), \
        int(f[7]), int(f[8]), int(f[9]), int(f[10]), int(f[11])
    if mq < MINQ or al < MINALN:
        continue
    side = q.rsplit("|", 1)[1]
    un = ql - qe if side == "L" else qs  # unaligned flank bases on the site side
    if un > SLACK:
        continue
    if side == "L":
        p = te + un if st == "+" else ts - un
    else:
        p = ts - un if st == "+" else te + un
    if q not in best or nm > best[q][0]:
        best[q] = (nm, t, max(0, min(tl, p)), st)
print("site\tstatus\tdist\tcontigL\tpointL\tstrandL\tcontigR\tpointR\tstrandR\tendL\tendR")
for l in list(open(sl))[1:]:
    s = l.split("\t")[0].strip()
    L, R = best.get(f"{s}|L"), best.get(f"{s}|R")
    atend = lambda x: int(min(x[2], ln[x[1]] - x[2]) <= END) if x else ""
    if L and R:
        d = R[2] - L[2] if L[3] == "+" else L[2] - R[2]
        if L[1] == R[1] and L[3] == R[3] and -SLACK <= d <= MAXDIST:
            st, d = "spanned", max(0, d)
        else:
            st, d = ("split_end" if atend(L) or atend(R) else "split_other"), ""
    else:
        st, d = ("left_only" if L else "right_only" if R else "none"), ""
    v = lambda x, i: x[i] if x else ""
    print(f"{s}\t{st}\t{d}\t{v(L, 1)}\t{v(L, 2)}\t{v(L, 3)}\t{v(R, 1)}\t{v(R, 2)}\t{v(R, 3)}\t{atend(L)}\t{atend(R)}")
