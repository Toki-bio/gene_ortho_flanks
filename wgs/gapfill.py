#!/usr/bin/env python3
"""Fill the N runs of a scaffolded assembly from a fresh haplotype-resolved assembly.

For every N run (>= MIN_N) of GENOME, the FLANK bp on each side (stopping at the next N run)
are aligned (minimap2 -c -x asm5) to each haplotype assembly separately, so MAPQ means
uniqueness within one haplotype. A contig fills the gap when it carries both flanks:
  - each flank alignment >= min(MIN_ALN, half the flank), MAPQ >= MIN_MAPQ;
  - same strand; left flank before right flank along the gap direction;
  - each alignment reaches within END_SLACK bp of the flank's gap-side end (no extrapolation
    beyond that: the contig position of the gap edge is the alignment end plus the unaligned
    flank bases, at most END_SLACK).
Per haplotype the best-supported pair (sum of matches) is kept.

Per haplotype status:
  spanned       fill = contig bases between the flanks (0 = flanks adjacent)
  split         both flanks found, on different contigs (or wrong order/strand)
  left_only / right_only / none
Gap class: filled (a haplotype has fill > JOIN_MAX), joined (all spanning haplotypes
<= JOIN_MAX), split, unspanned.

Writes OUT/gaps.tsv and OUT/fills.fa (>gapid|hapN|contig:start-end|strand, fill in GENOME
orientation, with FLANK_OUT bp of contig on each side for read checks).

Usage: gapfill.py GENOME.fa OUT_DIR LABEL=HAP.fa [LABEL=HAP.fa ...] [--threads N]
"""
import os
import re
import subprocess
import sys

MIN_N = 10
FLANK = 10000
MIN_ALN = 2000
MIN_MAPQ = 5
END_SLACK = 500
JOIN_MAX = 200
FLANK_OUT = 1000
comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")

args = [a for a in sys.argv[1:] if not a.startswith("--threads")]
threads = "16"
for i, a in enumerate(sys.argv):
    if a == "--threads":
        threads = sys.argv[i + 1]
        args.remove(threads)
genome, out = args[0], args[1]
haps = [a.split("=", 1) for a in args[2:]]
os.makedirs(out, exist_ok=True)


def read_fa(path):
    seqs, name, buf = {}, None, []
    for l in open(path):
        if l.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = l[1:].split()[0], []
        else:
            buf.append(l.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


g = read_fa(genome)
gaps = []                                   # (id, chrom, start, end, lflank(s,e), rflank(s,e))
for c, s in g.items():
    runs = [(m.start(), m.end()) for m in re.finditer("[Nn]{%d,}" % MIN_N, s)]
    for i, (a, b) in enumerate(runs):
        ls = max(runs[i - 1][1] if i else 0, a - FLANK)
        re_ = min(runs[i + 1][0] if i + 1 < len(runs) else len(s), b + FLANK)
        gaps.append((f"gap{len(gaps) + 1:04d}", c, a, b, (ls, a), (b, re_)))
fl = os.path.join(out, "flanks.fa")
with open(fl, "w") as h:
    for gid, c, a, b, (ls, le), (rs, re_) in gaps:
        if le - ls >= 1:
            h.write(f">{gid}|L\n{g[c][ls:le]}\n")
        if re_ - rs >= 1:
            h.write(f">{gid}|R\n{g[c][rs:re_]}\n")
print(f"{len(gaps)} N runs >= {MIN_N}", file=sys.stderr)


def hits(hapfa):
    o = subprocess.run(["minimap2", "-c", "-x", "asm5", "-t", threads, "--secondary=no", hapfa, fl],
                       capture_output=True, text=True, check=True).stdout
    h = {}
    for l in o.splitlines():
        f = l.split("\t")
        gid, side = f[0].split("|")
        h.setdefault(gid, {"L": [], "R": []})[side].append(dict(
            ql=int(f[1]), qs=int(f[2]), qe=int(f[3]), st=f[4], t=f[5], tl=int(f[6]), ts=int(f[7]),
            te=int(f[8]), match=int(f[9]), alen=int(f[10]), mapq=int(f[11])))
    return h


def ok(x):
    return x["qe"] - x["qs"] >= min(MIN_ALN, x["ql"] // 2) and x["mapq"] >= MIN_MAPQ


def edge(x, side):
    """Contig position of the gap edge of this flank and the unaligned slack; strand-aware.
    Left flank: its gap-side end is its 3' end (query end); right flank: its 5' end."""
    if side == "L":
        slack = x["ql"] - x["qe"]
        return (x["te"] + slack if x["st"] == "+" else x["ts"] - slack), slack
    slack = x["qs"]
    return (x["ts"] - slack if x["st"] == "+" else x["te"] + slack), slack


hapseq = {}
res = {}
for label, path in haps:
    hapseq[label] = read_fa(path)
    hh = hits(path)
    for gid, *_ in gaps:
        L = [x for x in hh.get(gid, {}).get("L", []) if ok(x)]
        R = [x for x in hh.get(gid, {}).get("R", []) if ok(x)]
        best = None
        for a in L:
            pa, sa = edge(a, "L")
            if sa > END_SLACK:
                continue
            for b in R:
                pb, sb = edge(b, "R")
                if b["t"] != a["t"] or b["st"] != a["st"] or sb > END_SLACK:
                    continue
                lo, hi = (pa, pb) if a["st"] == "+" else (pb, pa)
                if hi < lo - 50:              # flanks overlap by more than slack: wrong order
                    continue
                sup = a["match"] + b["match"]
                if best is None or sup > best[0]:
                    ident = (a["match"] + b["match"]) / (a["alen"] + b["alen"])
                    best = (sup, a["t"], max(0, lo), max(lo, hi), a["st"], ident)
        if best:
            res[(gid, label)] = ("spanned",) + best[1:]
        elif L and R:
            res[(gid, label)] = ("split", f"{L[0]['t']}/{R[0]['t']}")
        elif L:
            res[(gid, label)] = ("left_only", L[0]["t"])
        elif R:
            res[(gid, label)] = ("right_only", R[0]["t"])
        else:
            res[(gid, label)] = ("none",)

cols = ["gap", "chrom", "start", "end", "n_len", "class"]
for label, _ in haps:
    cols += [f"{label}_status", f"{label}_contig", f"{label}_fill", f"{label}_flank_ident"]
with open(os.path.join(out, "gaps.tsv"), "w") as h, open(os.path.join(out, "fills.fa"), "w") as fa:
    h.write("\t".join(cols) + "\n")
    for gid, c, a, b, *_ in gaps:
        row = [gid, c, a, b, b - a]
        per = []
        fills = []
        for label, _ in haps:
            r = res[(gid, label)]
            if r[0] == "spanned":
                _, t, lo, hi, st, ident = r
                fills.append(hi - lo)
                per += ["spanned", f"{t}:{lo}-{hi}", hi - lo, f"{ident:.4f}"]
                cs = hapseq[label][t]
                piece = cs[max(0, lo - FLANK_OUT):hi + FLANK_OUT]
                if st == "-":
                    piece = piece.translate(comp)[::-1]
                fa.write(f">{gid}|{label}|{t}:{lo}-{hi}|{st}|fill={hi - lo}\n{piece}\n")
            else:
                per += [r[0], r[1] if len(r) > 1 else "", "", ""]
        if fills:
            cls = "filled" if max(fills) > JOIN_MAX else "joined"
        elif any(res[(gid, l)][0] == "split" for l, _ in haps):
            cls = "split"
        else:
            cls = "unspanned"
        h.write("\t".join(map(str, row + [cls] + per)) + "\n")
print("done", file=sys.stderr)
