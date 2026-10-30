#!/usr/bin/env python3
"""Alignments for genes called empty_site: reference locus (gene + PAD each side)
vs the query site between the placed flanks (+ PAD), plus a row with only the
reference CDS bases (union of the gene's transcripts) to show where the gene was.

Usage: make_aln.py OUT_DIR(of gof.py) REF_BED ISOFORMS REF_2BIT QUERY_2BIT QUERY_ANNOT_BED GENE [GENE ...]
Writes OUT_DIR/aln/<gene>.aln.fa and <gene>.bed (ViewAlign ?bed= track, coordinates on
the named rows: flanks, reference transcripts as BED12 = exons/introns/CDS, query
TOGA2 projections of the gene).
"""
import os
import subprocess
import sys

PAD = 10000             # shows both flank scales of gof.py (1 kb, 10 kb)
FLANKS = (1000, 10000)
MAX_TX = 4              # distinct exon structures drawn per gene (ViewAlign has 6 lanes)
out, ref_bed, isoforms, r2bit, q2bit, q_annot = sys.argv[1:7]
genes = sys.argv[7:]
os.makedirs(os.path.join(out, "aln"), exist_ok=True)

comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def fetch(tb, c, s, e):
    r = subprocess.run(["twoBitToFa", f"{tb}:{c}:{s}-{e}", "stdout"], capture_output=True, text=True, check=True)
    return "".join(r.stdout.splitlines()[1:])


def bed_on_row(f, row, origin, rowlen, rc, rgb=None):
    """BED12 line f (genome coords) -> line on alignment row `row` (row coords)."""
    s, e = int(f[1]) - origin, int(f[2]) - origin
    if s < 0 or e > rowlen:
        return None
    blocks = [(s + int(b), s + int(b) + int(z)) for z, b in
              zip(f[10].rstrip(",").split(","), f[11].rstrip(",").split(","))]
    ts, te = int(f[6]) - origin, int(f[7]) - origin
    strand = f[5]
    if rc:
        def flip(x):
            return rowlen - x
        s, e = flip(e), flip(s)
        blocks = sorted((flip(b), flip(a)) for a, b in blocks)
        ts, te = flip(te), flip(ts)
        strand = {"+": "-", "-": "+"}.get(strand, strand)
    ts, te = max(s, ts), min(e, te)
    return "\t".join(map(str, [row, s, e, f[3], 0, strand, ts, te, rgb or f[8], len(blocks),
                               ",".join(str(b - a) for a, b in blocks) + ",",
                               ",".join(str(a - s) for a, b in blocks) + ","]))


tx2gene = dict(reversed(l.rstrip("\n").split("\t")[:2]) for l in open(isoforms))
cds = {g: set() for g in genes}
ref_tx = {g: [] for g in genes}
for line in open(ref_bed):
    f = line.rstrip("\n").split("\t")
    g = tx2gene.get(f[3])
    if g not in cds:
        continue
    key = f[1:3] + f[6:8] + f[10:12]
    if all(key != t[1:3] + t[6:8] + t[10:12] for t in ref_tx[g]):
        ref_tx[g].append(f)
    s0, ts, te = int(f[1]), int(f[6]), int(f[7])
    for size, st in zip(f[10].rstrip(",").split(","), f[11].rstrip(",").split(",")):
        a, b = s0 + int(st), s0 + int(st) + int(size)
        cds[g].update(range(max(a, ts), min(b, te)))

q_proj = {g: [] for g in genes}
for line in open(q_annot):
    f = line.rstrip("\n").split("\t")
    g = tx2gene.get("#".join(f[3].split("#")[:2]))
    if g in q_proj:
        q_proj[g].append(f)

rows = {}
hdr = None
for line in open(os.path.join(out, "verdicts.tsv")):
    f = line.rstrip("\n").split("\t")
    if hdr is None:
        hdr = f
        continue
    rows[f[0]] = dict(zip(hdr, f))

for g in genes:
    r = rows[g]
    c, rest = r["ref_locus"].split(":")
    s, e = map(int, rest[:-1].split("-"))
    q = r["f1000_query"] or r["f10000_query"]
    qc, qrest = q.split(":")
    qs, qe = map(int, qrest[:-1].split("-"))
    qstrand = qrest[-1]
    rs, re_ = s - PAD, e + PAD
    ref = fetch(r2bit, c, rs, re_).upper()
    qry = fetch(q2bit, qc, qs - PAD, qe + PAD).upper()
    if qstrand == "-":
        qry = qry.translate(comp)[::-1]
    raw = os.path.join(out, "aln", g + ".in.fa")
    rn = f"betulina_{c}:{rs}-{re_}_{g}_locus"
    qn = f"trizona_{qc}:{qs - PAD}-{qe + PAD}{'_rc' if qstrand == '-' else ''}_site"
    with open(raw, "w") as h:
        h.write(f">{rn}\n{ref}\n>{qn}\n{qry}\n")
    # genafpair + ep 0: one long gap for the missing gene instead of scattered gaps
    aln = subprocess.run(["mafft", "--genafpair", "--ep", "0", "--maxiterate", "1000", "--quiet", raw],
                         capture_output=True, text=True, check=True).stdout
    seqs, name = {}, None
    for l in aln.splitlines():
        if l.startswith(">"):
            name = l[1:]
            seqs[name] = []
        else:
            seqs[name].append(l.strip().upper())
    a_ref, a_qry = "".join(seqs[rn]), "".join(seqs[qn])
    cds_row, pos = [], rs
    for ch in a_ref:
        if ch == "-":
            cds_row.append("-")
            continue
        cds_row.append(ch if pos in cds[g] else "-")
        pos += 1
    dst = os.path.join(out, "aln", g + ".aln.fa")
    with open(dst, "w") as h:
        h.write(f">{rn}\n{a_ref}\n>betulina_{g}_CDS_only\n{''.join(cds_row)}\n>{qn}\n{a_qry}\n")
    os.remove(raw)

    # annotation track on the named rows
    bed = [f"track name={g}_flanks_exons"]
    rlen = re_ - rs
    for fl in FLANKS:
        for side, a, b in (("L", PAD - fl, PAD), ("R", rlen - PAD, rlen - PAD + fl)):
            bed.append(f"{rn}\t{a}\t{b}\t{side}_flank_{fl // 1000}kb\t0\t.\t{a}\t{a}\t150,150,150")
    for f in ref_tx[g][:MAX_TX]:
        line = bed_on_row(f, rn, rs, rlen, False, "0,0,150")
        if line:
            bed.append(line)
    for f in q_proj[g]:
        if f[0] == qc:
            line = bed_on_row(f, qn, qs - PAD, len(qry), qstrand == "-")
            if line:
                bed.append(line)
    with open(os.path.join(out, "aln", g + ".bed"), "w") as h:
        h.write("\n".join(bed) + "\n")
    print(dst, len(a_ref), "columns,", len(bed) - 1, "BED features")
