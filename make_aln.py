#!/usr/bin/env python3
"""Alignments for genes called empty_site: reference locus (gene + PAD each side)
vs the query site between the placed flanks (+ PAD), plus a row with only the
reference CDS bases (union of the gene's transcripts) to show where the gene was.

Usage: make_aln.py OUT_DIR(of gof.py) REF_BED ISOFORMS REF_2BIT QUERY_2BIT GENE [GENE ...]
Writes OUT_DIR/aln/<gene>.aln.fa
"""
import os
import subprocess
import sys

PAD = 2000
out, ref_bed, isoforms, r2bit, q2bit = sys.argv[1:6]
genes = sys.argv[6:]
os.makedirs(os.path.join(out, "aln"), exist_ok=True)

comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def fetch(tb, c, s, e):
    r = subprocess.run(["twoBitToFa", f"{tb}:{c}:{s}-{e}", "stdout"], capture_output=True, text=True, check=True)
    return "".join(r.stdout.splitlines()[1:])


tx2gene = dict(reversed(l.rstrip("\n").split("\t")[:2]) for l in open(isoforms))
cds = {g: set() for g in genes}
for line in open(ref_bed):
    f = line.rstrip("\n").split("\t")
    g = tx2gene.get(f[3])
    if g not in cds:
        continue
    s0, ts, te = int(f[1]), int(f[6]), int(f[7])
    for size, st in zip(f[10].rstrip(",").split(","), f[11].rstrip(",").split(",")):
        a, b = s0 + int(st), s0 + int(st) + int(size)
        cds[g].update(range(max(a, ts), min(b, te)))

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
    print(dst, len(a_ref), "columns")
