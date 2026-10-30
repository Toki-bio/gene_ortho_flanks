#!/usr/bin/env python3
"""Evaluate per-region HiFi rebuilds.

For each region (regions.tsv) with contigs in ASM/<id>/contigs.fa:
  - contigs are aligned (minimap2 -c asm5) to the padded query bait <id>|q, whose
    region interval is [rlo, rhi) (gap = N run, deletion = empty site, control = window);
  - a contig spans the region when its alignment(s) cover both rlo - M and rhi + M
    (M = 500, so it does not matter where the aligner puts an indel at the edge);
    bait positions are converted to contig positions through the CIGAR, because one
    alignment can carry a large insertion (toy: a 26.7 kb insertion in one block);
  - rebuilt = contig bases between the two points minus 2M; it is aligned to the
    reference bait <id>|r: r_cov = reference bases covered, r_ident = identity.
Status: spanned / one_side / no_span / no_contig. Controls: rebuilt vs the assembly.

Usage: evaluate.py regions.tsv bait.fa ASM_DIR PAD
"""
import os
import re
import subprocess
import sys
import tempfile

regions_f, bait_f, asm, PAD = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
M = 500                  # measure this far outside each region edge
MAX_END_DIST = 2000      # an alignment may stop this short of a measuring point (extrapolated)


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


def paf(target_seq, query_seqs, preset="asm5"):
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as d:
        t, q = os.path.join(d, "t.fa"), os.path.join(d, "q.fa")
        open(t, "w").write(f">t\n{target_seq}\n")
        open(q, "w").write("".join(f">{n}\n{s}\n" for n, s in query_seqs.items()))
        out = subprocess.run(["minimap2", "-c", "-x", preset, "--secondary=no", t, q],
                             capture_output=True, text=True, check=True).stdout
    rows = []
    for l in out.splitlines():
        f = l.split("\t")
        cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
        rows.append(dict(q=f[0], ql=int(f[1]), qs=int(f[2]), qe=int(f[3]), st=f[4], ts=int(f[7]), te=int(f[8]),
                         match=int(f[9]), alen=int(f[10]), cigar=re.findall(r"(\d+)([MIDNSHP=X])", cg)))
    return rows


def t2q(a, t):
    """Contig (PAF query) position aligned to bait (PAF target) position t; outside the
    alignment, extrapolated from its nearer end."""
    if t <= a["ts"]:
        off = -(a["ts"] - t)
    else:
        tpos, off = a["ts"], 0
        for n, op in a["cigar"]:
            n = int(n)
            if op in "M=X":
                if t < tpos + n:
                    off += t - tpos
                    break
                tpos += n
                off += n
            elif op == "I":
                off += n
            elif op in "DN":
                if t < tpos + n:
                    break
                tpos += n
        else:
            off += t - tpos
    return a["qs"] + off if a["st"] == "+" else a["qe"] - off


def reaches(a, t):
    return a["ts"] - MAX_END_DIST <= t <= a["te"] + MAX_END_DIST


bait = read_fa(bait_f)
print("id\tkind\tgenes\treads\tcontigs\tstatus\tquery_len\trebuilt_len\tr_len\tr_cov\tr_ident\tnote")
for l in open(regions_f):
    f = l.rstrip("\n").split("\t")
    if f[0] == "id":
        continue
    rid, kind, qc, qs, qe, rc, rs, re_, genes = f[0], f[1], f[2], int(f[3]), int(f[4]), f[5], int(f[6]), int(f[7]), f[8]
    off = qs - max(0, qs - PAD)           # region start inside the query bait
    rlo, rhi = off, off + (qe - qs)
    d = os.path.join(asm, rid)
    names = os.path.join(d, "names.txt")
    nreads = sum(1 for _ in open(names)) if os.path.exists(names) else 0
    cf = os.path.join(d, "contigs.fa")
    contigs = read_fa(cf) if os.path.exists(cf) else {}
    base = [rid, kind, genes, str(nreads), str(len(contigs))]
    if not contigs:
        print("\t".join(base + ["no_contig", str(qe - qs), "", str(re_ - rs), "", "", ""]))
        continue
    aln = paf(bait[f"{rid}|q"], contigs)
    best = None
    for name in contigs:
        a = [x for x in aln if x["q"] == name]
        for L in a:
            if not reaches(L, rlo - M):
                continue
            for R in a:
                if R["st"] != L["st"] or not reaches(R, rhi + M) or R["te"] < L["te"]:
                    continue
                c1, c2 = t2q(L, rlo - M), t2q(R, rhi + M)
                lo, hi = sorted((c1, c2))
                if lo < 0 or hi > len(contigs[name]):
                    continue
                if best is None or (hi - lo) > (best[2] - best[1]):
                    best = (name, lo, hi)
    if best is None:
        touched = any(x["te"] >= rlo - MAX_END_DIST and x["ts"] <= rhi + MAX_END_DIST for x in aln)
        print("\t".join(base + ["one_side" if touched else "no_span", str(qe - qs), "", str(re_ - rs), "", "", ""]))
        continue
    name, lo, hi = best
    middle = contigs[name][lo + M:hi - M] if hi - lo > 2 * M else ""
    note = r_cov = r_ident = ""
    if len(middle) >= 50:
        ra = paf(bait[f"{rid}|r"], {"middle": middle}, "asm20")
        cov = set()
        m = al = 0
        for x in ra:
            cov.update(range(x["ts"], x["te"]))
            m += x["match"]
            al += x["alen"]
        r_cov = f"{len(cov) / max(1, re_ - rs):.3f}"
        r_ident = f"{m / al:.4f}" if al else ""
    if kind == "control":
        ca = paf(bait[f"{rid}|q"][rlo:rhi], {"middle": middle}, "asm5") if middle else []
        m = sum(x["match"] for x in ca)
        al = sum(x["alen"] for x in ca)
        note = f"vs_assembly_ident={m / al:.4f}" if al else "vs_assembly_unaligned"
    print("\t".join(base + ["spanned", str(qe - qs), str(hi - lo - 2 * M), str(re_ - rs), r_cov, r_ident, note]))
