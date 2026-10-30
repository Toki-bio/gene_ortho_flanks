#!/usr/bin/env python3
"""ViewAlign alignments of rebuilt regions: reference (betulina) stretch + PAD each side,
the query (trizona) assembly around the region (N run / empty site), and the rebuilt
haplotype contigs (spans/<rid>__<contig>.fa from evaluate.py SPAN_DIR), all in the query
assembly's orientation; MAFFT (--auto).

BED (?bed=): reference transcripts of the region's genes (BED12, exons/CDS) on the reference
row, the region (N run / empty site) on the assembly row, the rebuilt middle on each contig row.

Usage: aln_rebuild.py regions.tsv bait.fa SPAN_DIR REF_2BIT REF_BED ISOFORMS OUT_DIR RID [RID ...]
"""
import glob
import os
import subprocess
import sys
import tempfile

regions_f, bait_f, sdir, r2bit, ref_bed, isoforms, out = sys.argv[1:8]
rids = sys.argv[8:]
PAD = 10000
os.makedirs(out, exist_ok=True)
comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(comp)[::-1]


def read_fa(path):
    seqs, name = {}, None
    for l in open(path):
        l = l.strip()
        if l.startswith(">"):
            name = l[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(l)
    return {k: "".join(v) for k, v in seqs.items()}


def fetch(tb, c, s, e):
    r = subprocess.run(["twoBitToFa", f"{tb}:{c}:{s}-{e}", "stdout"], capture_output=True, text=True, check=True)
    return "".join(r.stdout.splitlines()[1:])


def strand_vs(target, query, preset):
    """+ or - : strand of query's longest alignment on target."""
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as d:
        t, q = os.path.join(d, "t.fa"), os.path.join(d, "q.fa")
        open(t, "w").write(f">t\n{target}\n")
        open(q, "w").write(f">q\n{query}\n")
        o = subprocess.run(["minimap2", "-x", preset, t, q], capture_output=True, text=True, check=True).stdout
    best = max((l.split("\t") for l in o.splitlines()), key=lambda f: int(f[10]), default=None)
    return best[4] if best else "+"


def bed12_on_row(f, row, origin, rowlen, flip, rgb):
    s, e = int(f[1]) - origin, int(f[2]) - origin
    if s < 0 or e > rowlen:
        return None
    blocks = [(s + int(b), s + int(b) + int(z)) for z, b in
              zip(f[10].rstrip(",").split(","), f[11].rstrip(",").split(","))]
    ts, te, strand = int(f[6]) - origin, int(f[7]) - origin, f[5]
    if flip:
        s, e = rowlen - e, rowlen - s
        blocks = sorted((rowlen - b, rowlen - a) for a, b in blocks)
        ts, te = rowlen - te, rowlen - ts
        strand = {"+": "-", "-": "+"}[strand]
    ts, te = max(s, ts), min(e, te)
    return "\t".join(map(str, [row, s, e, f[3], 0, strand, ts, te, rgb, len(blocks),
                               ",".join(str(b - a) for a, b in blocks) + ",",
                               ",".join(str(a - s) for a, b in blocks) + ","]))


tx2gene = dict(reversed(l.rstrip("\n").split("\t")[:2]) for l in open(isoforms))
bait = read_fa(bait_f)
regions = {}
for l in open(regions_f):
    f = l.rstrip("\n").split("\t")
    if f[0] != "id":
        regions[f[0]] = f
for rid in rids:
    _, kind, qc, qs, qe, rc_, rs, re_, genes = regions[rid]
    qs, qe, rs, re_ = int(qs), int(qe), int(rs), int(re_)
    q = bait[f"{rid}|q"]                                  # query region +- PAD
    qoff = qs - max(0, qs - PAD)
    r0 = max(0, rs - PAD)
    r = fetch(r2bit, rc_, r0, re_ + PAD)
    rflip = strand_vs(q, r, "asm20") == "-"
    if rflip:
        r = rc(r)
    rows = [(f"betulina_{rc_}_{r0}", r), (f"trizona_assembly_{qc}_{qs - qoff}", q)]
    bed = []
    gset = set(genes.split(","))
    seen = set()
    for line in open(ref_bed):
        f = line.rstrip("\n").split("\t")
        if f[0] != rc_ or tx2gene.get(f[3]) not in gset:
            continue
        key = tuple(f[1:3] + f[6:8] + f[10:12])
        if key in seen or len(seen) >= 4:
            continue
        seen.add(key)
        b = bed12_on_row(f, rows[0][0], r0, len(r), rflip, "0,0,150")
        if b:
            bed.append(b)
    rlo, rhi = rs - r0, re_ - r0
    if rflip:
        rlo, rhi = len(r) - rhi, len(r) - rlo
    bed.append(f"{rows[0][0]}\t{rlo}\t{rhi}\treference stretch\t0\t.\t{rlo}\t{rlo}\t150,150,150")
    label = "N run" if kind == "gap" else "query site" if kind == "deletion" else "window"
    bed.append(f"{rows[1][0]}\t{qoff}\t{qoff + qe - qs}\t{label} ({qe - qs} bp)\t0\t.\t{qoff}\t{qoff}\t200,0,0")
    for fpath in sorted(f for f in glob.glob(os.path.join(sdir, f"{rid}__*.fa")) if not f.endswith(".mid.fa")):
        hdr, seq = open(fpath).read().split("\n")[:2]
        name = hdr[1:].split()[0].split("__")[1]
        a, b = map(int, hdr.split("middle=")[1].split("-"))
        if strand_vs(q, seq, "asm5") == "-":
            seq = rc(seq)
            a, b = len(seq) - b, len(seq) - a
        rows.append((f"rebuilt_{name}", seq))
        bed.append(f"rebuilt_{name}\t{a}\t{max(b, a + 1)}\trebuilt ({b - a} bp)\t0\t.\t{a}\t{a}\t0,140,0")
    raw = os.path.join(out, f"{rid}.raw.fa")
    open(raw, "w").write("".join(f">{n}\n{s}\n" for n, s in rows))
    aln = subprocess.run(["mafft", "--auto", "--thread", "16", "--quiet", raw], capture_output=True, text=True, check=True).stdout
    open(os.path.join(out, f"{rid}.aln.fa"), "w").write(aln)
    open(os.path.join(out, f"{rid}.bed"), "w").write("\n".join(bed) + "\n")
    os.remove(raw)
    print(rid, kind, genes, [(n, len(s)) for n, s in rows], "rflip" if rflip else "")
