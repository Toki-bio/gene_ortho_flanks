#!/usr/bin/env python3
"""What do reads hold right after the gap-side end of each flank? (step 1 of the gap plan)

FILLS.fa: one entry per region, header `>NAME middle=A-B`: left flank = [0,A), piece = [A,B), right flank = [B,end)
(entries named `mid` are skipped). Reads are aligned (minimap2 map-hifi) to the three parts as separate targets, then
followed across each gap end in the region's forward direction:
  from the left flank's gap end (reads reaching it):   piece  = next aligned part starts the piece
                                                       direct = next aligned part starts the right flank
                                                       other  = >= MIN_CONT read bases follow, aligned to none of the parts
                                                       elsewhere = the next aligned part starts inside a part, not at its start
                                                       ends   = fewer than MIN_CONT read bases follow
  from the right flank's gap start, the same backwards (piece = previous part ends the piece).
Writes OUT/ends.tsv (counts per region, side, class), OUT/ends_reads.tsv (read, region, side, class, read position)
and OUT/other_tails.fa (the first TAIL bp after the gap end of every `other` read, in the region's direction).
Usage: classify_ends.py FILLS.fa READS.fq OUT_DIR [THREADS]"""
import os
import subprocess
import sys
from collections import Counter, defaultdict

EDGE = 200        # an alignment "reaches" a part end when it stops within EDGE bp of it
MIN_CONT = 1000   # read bases needed beyond the gap end to say what follows
JOIN_SLACK = 300  # next alignment must begin within this many read bases of the gap position
TAIL = 3000
comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
rc = lambda s: s.translate(comp)[::-1]

fills, reads, out = sys.argv[1:4]
threads = sys.argv[4] if len(sys.argv) > 4 else "8"
os.makedirs(out, exist_ok=True)


def read_fa(path):
    seqs, name, head, buf = {}, None, {}, []
    for l in open(path):
        if l.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = l[1:].split()[0], []
            head[name] = l[1:].strip()
        else:
            buf.append(l.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs, head


seqs, head = read_fa(fills)
tlen = {}
with open(os.path.join(out, "parts.fa"), "w") as fh:
    for n, s in seqs.items():
        if n == "mid":
            continue
        a, b = map(int, head[n].split("middle=")[1].split()[0].split("-"))
        for part, x in (("L", s[:a]), ("M", s[a:b]), ("R", s[b:])):
            if x:
                fh.write(f">{n}|{part}\n{x}\n")
                tlen[(n, part)] = len(x)

segs = defaultdict(list)          # (read, region) -> [(fs, fe, part, ts, te)] in the region's forward read frame
qlen = {}
p = subprocess.Popen(["minimap2", "-c", "-x", "map-hifi", "--secondary=no", "-t", threads,
                      os.path.join(out, "parts.fa"), reads], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
for l in p.stdout:
    f = l.split("\t")
    reg, part = f[5].rsplit("|", 1)
    q, ql, qs, qe, st, ts, te = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], int(f[7]), int(f[8])
    qlen[q] = ql
    fs, fe = (qs, qe) if st == "+" else (ql - qe, ql - qs)
    segs[(q, reg)].append((fs, fe, part, ts, te, st))
if p.wait():
    sys.exit("minimap2 failed")


def strand_ok(a, b):
    return a[5] == b[5]


def covered(sg, x, y):
    return any(s[0] < y and s[1] > x for s in sg)


cnt = Counter()
rows = []
need_tail = {}
for (q, reg), sg in segs.items():
    ql = qlen[q]
    # from the left flank's gap end, forwards
    for s in sg:
        if s[2] != "L" or s[4] < tlen[(reg, "L")] - EDGE:
            continue
        pos = s[1] + (tlen[(reg, "L")] - s[4])
        if ql - pos < MIN_CONT:
            cls = "ends"
        else:
            nxt = [t for t in sg if t is not s and strand_ok(s, t) and abs(t[0] - pos) <= JOIN_SLACK]
            if any(t[2] == "M" and t[3] <= EDGE for t in nxt):
                cls = "piece"
            elif any(t[2] == "R" and t[3] <= EDGE for t in nxt):
                cls = "direct"
            elif nxt:
                cls = "elsewhere:" + ",".join(f"{t[2]}{t[3]}" for t in nxt)
            elif not covered(sg, pos + 50, pos + MIN_CONT):
                cls = "other"
                need_tail[(q, reg, "L")] = (pos, s[5])
            else:
                cls = "elsewhere"
        cnt[(reg, "L", cls.split(":")[0])] += 1
        rows.append((q, reg, "L", cls, pos))
    # from the right flank's gap start, backwards
    for s in sg:
        if s[2] != "R" or s[3] > EDGE:
            continue
        pos = s[0] - s[3]
        if pos < MIN_CONT:
            cls = "ends"
        else:
            prv = [t for t in sg if t is not s and strand_ok(s, t) and abs(t[1] - pos) <= JOIN_SLACK]
            if any(t[2] == "M" and t[4] >= tlen[(reg, "M")] - EDGE for t in prv):
                cls = "piece"
            elif any(t[2] == "L" and t[4] >= tlen[(reg, "L")] - EDGE for t in prv):
                cls = "direct"
            elif prv:
                cls = "elsewhere:" + ",".join(f"{t[2]}{t[4]}" for t in prv)
            elif not covered(sg, pos - MIN_CONT, pos - 50):
                cls = "other"
                need_tail[(q, reg, "R")] = (pos, s[5])
            else:
                cls = "elsewhere"
        cnt[(reg, "R", cls.split(":")[0])] += 1
        rows.append((q, reg, "R", cls, pos))

# tails of `other` reads, in the region's forward direction (R side: the TAIL bp before the gap start)
want = {q for q, _, _ in need_tail}
with open(os.path.join(out, "other_tails.fa"), "w") as fh:
    with open(reads) as rf:
        while True:
            h = rf.readline()
            if not h:
                break
            s = rf.readline().strip()
            rf.readline(); rf.readline()
            q = h[1:].split()[0]
            if q not in want:
                continue
            for (qq, reg, side), (pos, st) in need_tail.items():
                if qq != q:
                    continue
                fwd = s if st == "+" else rc(s)
                t = fwd[pos:pos + TAIL] if side == "L" else fwd[max(0, pos - TAIL):pos]
                fh.write(f">{q}|{reg}|{side}\n{t}\n")

classes = ["piece", "direct", "other", "elsewhere", "ends"]
regs = sorted({r for r, _, _ in cnt})
with open(os.path.join(out, "ends.tsv"), "w") as fh:
    fh.write("region\tside\t" + "\t".join(classes) + "\n")
    for r in regs:
        for side in ("L", "R"):
            fh.write(f"{r}\t{side}\t" + "\t".join(str(cnt[(r, side, c)]) for c in classes) + "\n")
with open(os.path.join(out, "ends_reads.tsv"), "w") as fh:
    fh.write("read\tregion\tside\tclass\tread_pos\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
print(open(os.path.join(out, "ends.tsv")).read(), end="")
