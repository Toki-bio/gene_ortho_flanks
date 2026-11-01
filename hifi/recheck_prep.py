#!/usr/bin/env python3
"""Re-examination of the weaker gap labels (owner's challenges, 2026-11-01): build read-fishing windows.

For every gap of the selected causes (gap_causes.tsv column `cause`) write
  win.fa    the deposited hap1 window [start - W, end + W] (N run kept), named `gapNNNN|h1|chrom:a-b|gap=x-y`
            (x-y = N run in window coordinates); for gaps with a fill or join in another assembly also the alternative
            window (fill + W each side), named `gapNNNN|alt|SOURCE|contig:a-b|fill=x-y` (dep2 first, else asm25 h2, h1)
  cand.fa   25-mers every STEP bp of all windows (no N), to be counted over hap1 + hap2 (jellyfish) -> words step
  targets.tsv  gap, cause, chrom, start, end, alt source/contig/fill coordinates
Usage:
  recheck_prep.py prep HAP1.fa GAPS.tsv CAUSES.tsv CAUSE[,CAUSE..] EXTRA_GAPS(g,g|-) OUT [W] [STEP]
     (other assemblies: DEP2 / ASM_H1 / ASM_H2 FASTA paths from env DEP2, ASM_H1, ASM_H2; each needs a .fai)
  recheck_prep.py words OUT COUNTS.txt MAXCOUNT   (jellyfish dump -c) -> OUT/anchors.txt"""
import os
import sys

K = 25


def fai(path):
    d = {}
    for l in open(path + ".fai"):
        f = l.split("\t")
        d[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    return d


def fetch(path, idx, chrom, a, b):
    ln, off, lb, lby = idx[chrom]
    a, b = max(0, a), min(ln, b)
    with open(path, "rb") as fh:
        o1 = off + (a // lb) * lby + a % lb
        o2 = off + (b // lb) * lby + b % lb
        fh.seek(o1)
        s = fh.read(o2 - o1).decode().replace("\n", "").replace("\r", "")
    return s.upper(), a, b


def prep(argv):
    hap1, gaps, causes, sel, extra, out = argv[:6]
    W = int(argv[6]) if len(argv) > 6 else 40000
    STEP = int(argv[7]) if len(argv) > 7 else 15
    os.makedirs(out, exist_ok=True)
    sel = set(sel.split(","))
    extra = set() if extra == "-" else set(extra.split(","))
    cause = {}
    for l in list(open(causes))[1:]:
        f = l.rstrip("\n").split("\t")
        cause[f[0]] = f[-1]
    asm = {"dep2": os.environ.get("DEP2"), "asm_h2": os.environ.get("ASM_H2"), "asm_h1": os.environ.get("ASM_H1")}
    idx = {k: fai(v) for k, v in asm.items() if v}
    i1 = fai(hap1)
    rows = list(open(gaps))
    hd = rows[0].rstrip("\n").split("\t")
    win = open(os.path.join(out, "win.fa"), "w")
    cand = open(os.path.join(out, "cand.fa"), "w")
    tg = open(os.path.join(out, "targets.tsv"), "w")
    tg.write("gap\tcause\tchrom\tstart\tend\talt_source\talt_contig\talt_fill\n")
    nw = nc = 0

    def words(name, s):
        nonlocal nc
        for i in range(0, len(s) - K + 1, STEP):
            w = s[i:i + K]
            if "N" not in w:
                cand.write(f">{nc}\n{w}\n"); nc += 1

    for l in rows[1:]:
        f = dict(zip(hd, l.rstrip("\n").split("\t")))
        g = f["gap"]
        if cause.get(g) not in sel and g not in extra:
            continue
        c, s, e = f["chrom"], int(f["start"]), int(f["end"])
        seq, a, b = fetch(hap1, i1, c, s - W, e + W)
        name = f"{g}|h1|{c}:{a}-{b}|gap={s - a}-{e - a}"
        win.write(f">{name}\n{seq}\n"); words(name, seq); nw += 1
        alt = ("", "", "")
        for src in ("dep2", "asm_h2", "asm_h1"):
            st, ctg = f.get(f"{src}_status", ""), f.get(f"{src}_contig", "")
            if st == "spanned" and ":" in ctg and src in idx:
                cc, xy = ctg.rsplit(":", 1)
                x, y = map(int, xy.split("-"))
                aseq, aa, bb = fetch(asm[src], idx[src], cc, x - W, y + W)
                an = f"{g}|alt|{src}|{cc}:{aa}-{bb}|fill={x - aa}-{y - aa}"
                win.write(f">{an}\n{aseq}\n"); words(an, aseq); nw += 1
                alt = (src, ctg, str(y - x))
                break
        tg.write(f"{g}\t{cause.get(g)}\t{c}\t{s}\t{e}\t" + "\t".join(alt) + "\n")
    print(f"{nw} windows, {nc} candidate words", file=sys.stderr)


def words_step(argv):
    out, counts, mx = argv[0], argv[1], int(argv[2])
    keep = 0
    with open(os.path.join(out, "anchors.txt"), "w") as fh:
        for l in open(counts):
            w, n = l.split()
            if 1 <= int(n) <= mx:
                fh.write(w + "\n"); keep += 1
    print(f"{keep} anchor words", file=sys.stderr)


if __name__ == "__main__":
    {"prep": prep, "words": words_step}[sys.argv[1]](sys.argv[2:])
