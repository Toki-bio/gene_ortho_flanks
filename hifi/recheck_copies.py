#!/usr/bin/env python3
"""Where are the other copies of a collapsed / high-copy gap flank, and where do read tails that leave it go?

  prep:  recheck_copies.py prep WIN.fa CAUSES.tsv CAUSE[,CAUSE] OUT.fa [SEG]
         -> the SEG bp of hap1 flank next to each N run (`gap|L|chrom:pos` / `gap|R|chrom:pos`, pos = window start)
  table: recheck_copies.py table SEGS.paf TAILS.paf CAUSES.tsv CAUSE[,CAUSE] TAILS.fa > copies.tsv
         SEGS.paf: minimap2 -c -x asm20 -N 200 -p 0.1 of OUT.fa to hap1; TAILS.paf: the same for recheck_reads.py
         tails.fa (read sequence beyond the gap end, read name `gap|side|read`).
Per segment (index names may lack the .1 version suffix; both forms count as the own chromosome): loci (distinct chrom / 50 kb) with an alignment covering >= 50 % of the segment at >= 90 % identity, the
own locus excluded (same chrom within 100 kb of the gap); best identity of another copy; copies at >= 99 %.
Per gap side: tails mapped (>= 1 kb aligned, >= 95 % id): to the own locus (other side of the gap, i.e. the gap
is crossed), to another locus (named), or unmapped."""
import sys
from collections import Counter, defaultdict

if sys.argv[1] == "prep":
    win, causes, sel, out = sys.argv[2:6]
    SEG = int(sys.argv[6]) if len(sys.argv) > 6 else 3000
    sel = set(sel.split(","))
    want = {l.split("\t")[0] for l in open(causes) if l.rstrip("\n").split("\t")[-1] in sel}
    name, seq = None, {}
    for l in open(win):
        if l.startswith(">"):
            name = l[1:].strip(); seq[name] = []
        else:
            seq[name].append(l.strip())
    with open(out, "w") as fh:
        for n, s in seq.items():
            f = n.split("|")
            if f[1] != "h1" or f[0] not in want:
                continue
            s = "".join(s)
            x, y = map(int, f[3].split("=")[1].split("-"))
            c, ab = f[2].rsplit(":", 1)
            a = int(ab.split("-")[0])
            L = s[max(0, x - SEG):x]
            L = L[L.rfind("N") + 1:]
            R = s[y:y + SEG]
            R = R[:R.find("N")] if "N" in R else R
            if len(L) >= 500:
                fh.write(f">{f[0]}|L|{c}:{a + x}\n{L}\n")
            if len(R) >= 500:
                fh.write(f">{f[0]}|R|{c}:{a + y}\n{R}\n")
    sys.exit()

segpaf, tailpaf, causes, sel, tailfa = sys.argv[2:7]
alltails = defaultdict(int)
for l in open(tailfa):
    if l.startswith(">"):
        g, side, r = l[1:].strip().split("|", 2)
        alltails[(g, side)] += 1
sel = set(sel.split(","))
cause = {}
for l in list(open(causes))[1:]:
    f = l.rstrip("\n").split("\t")
    cause[f[0]] = (f[-1], f[1])
seg = defaultdict(list)
qlen = {}
for l in open(segpaf):
    f = l.split("\t")
    q, ql, qs, qe, t, ts = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[7])
    qlen[q] = ql
    idn = int(f[9]) / int(f[10])
    if (qe - qs) >= 0.5 * ql and idn >= 0.90:
        seg[q].append((t, ts, idn))
tails = defaultdict(list)
for l in open(tailpaf):
    f = l.split("\t")
    idn = int(f[9]) / int(f[10])
    if int(f[3]) - int(f[2]) >= 1000 and idn >= 0.95:
        tails[f[0]].append((f[5], int(f[7]), idn, int(f[11])))
print("gap\tside\tchrom\tseg_len\tother_loci\tloci_ge99\tbest_other_id\texamples\ttails\ttails_own_locus\ttails_other\ttails_unmapped\ttail_targets")
tail_by = defaultdict(lambda: defaultdict(list))
for q in tails:
    g, side, r = q.split("|", 2)
    tail_by[(g, side)][r] = tails[q]
segs = sorted(set(list(qlen)) | {q for q in seg})
for q in segs:
    g, side, pos = q.split("|")
    if cause.get(g, ("",))[0] not in sel:
        continue
    c, p = pos.rsplit(":", 1)
    p = int(p)
    loci = {}
    for t, ts, idn in seg[q]:
        if t in (c, c.rsplit(".", 1)[0]) and abs(ts - p) < 100000:
            continue
        k = (t, ts // 50000)
        loci[k] = max(loci.get(k, 0), idn)
    ex = ",".join(f"{t}:{k * 50}kb:{v:.3f}" for (t, k), v in sorted(loci.items(), key=lambda kv: -kv[1])[:3])
    tb = tail_by.get((g, side), {})
    own = oth = 0
    tg = Counter()
    for r, hits in tb.items():
        best = max(hits, key=lambda h: h[2])
        if best[0] in (c, c.rsplit(".", 1)[0]) and abs(best[1] - p) < 200000:
            own += 1
        else:
            oth += 1; tg[f"{best[0]}:{best[1] // 100000 / 10:.1f}Mb"] += 1
    print(f"{g}\t{side}\t{c}\t{qlen.get(q, 0)}\t{len(loci)}\t{sum(v >= 0.99 for v in loci.values())}\t"
          f"{max(loci.values(), default=0):.3f}\t{ex}\t{alltails[(g, side)]}\t{own}\t{oth}\t{alltails[(g, side)] - own - oth}\t{','.join(f'{k}x{v}' for k, v in tg.most_common(3))}")
