#!/usr/bin/env python3
"""gene_ortho_flanks: flank and neighbour-gene evidence for TOGA2 gene calls.

For each reference gene: map its flanks (FLANKS, e.g. 1 kb and 10 kb next to the
gene) to the query and place its nearest intact one2one neighbour genes, then
call per scale: empty_site / present / size_changed / assembly_gap /
assembly_edge / one_flank / synteny_broken / no_anchor.

Usage: gof.py --toga TOGA_OUT --ref-bed BED12 --isoforms TSV --ref-2bit R.2bit
              --query-2bit Q.2bit --query-index Q.mmi --genes GENES.txt --out DIR
"""
import argparse
import bisect
import collections
import os
import re
import struct
import subprocess
import sys

FLANKS = (1000, 10000)
MIN_MAPQ = 20
MIN_QCOV = 0.5          # aligned part of a flank ...
MIN_ALIGNED = 300       # ... or at least this many bp (a 1 kb flank often aligns only next to the gene)
EMPTY_FRAC = 0.2        # query gap <= 20 % of the gene span -> empty site
PRESENT_LO, PRESENT_HI = 0.5, 2.0
INTACT = ("FI", "I")
NONUNIQUE_S2S1 = 0.8    # flank whose best secondary chain scores >= 0.8 of the primary is not unique
                        # (minimap2's own -p default; b2t: 4 % of FI-control vs 93 % of IG/TR flanks)
EDGE_FLANK = 1000       # edge scale: flank length at each edge of the unaligned stretch
MAX_EDGE_RUN = 2_000_000
MIN_N_RUN = 10          # N run between anchors that makes the interval an assembly gap
MAX_SLACK = 100         # empty site needs aligned flank within this of the gene (or 20 % of its length)


def sh(cmd):
    r = subprocess.run(cmd, shell=True, executable="/bin/bash")
    if r.returncode:
        sys.exit("failed: " + cmd)


def read_genes(ref_bed, isoforms):
    tx2gene = {}
    with open(isoforms) as h:
        for line in h:
            g, t = line.rstrip("\n").split("\t")[:2]
            tx2gene[t] = g
    span = {}
    with open(ref_bed) as h:
        for line in h:
            f = line.rstrip("\n").split("\t")
            g = tx2gene.get(f[3])
            if g is None:
                continue
            c, s, e, st = f[0], int(f[1]), int(f[2]), f[5]
            if g in span:
                c0, s0, e0, st0 = span[g]
                if c0 != c:
                    continue
                span[g] = (c, min(s, s0), max(e, e0), st0)
            else:
                span[g] = (c, s, e, st)
    return span


def read_toga(toga):
    status = {}
    with open(os.path.join(toga, "loss_summary.tsv")) as h:
        for line in h:
            f = line.rstrip("\n").split("\t")
            if f[0] == "GENE":
                status[f[1]] = f[2]
    proj = {}
    with open(os.path.join(toga, "query_annotation.bed")) as h:
        for line in h:
            f = line.rstrip("\n").split("\t")
            proj[f[3]] = (f[0], int(f[1]), int(f[2]), f[5])
    o2o = {}
    with open(os.path.join(toga, "orthology_classification.tsv")) as h:
        next(h)
        for line in h:
            tg, tt, qg, qt, cls = line.rstrip("\n").split("\t")
            if cls == "one2one" and qt in proj and tg not in o2o:
                o2o[tg] = (tt, proj[qt])
    # all query projections per reference gene (any status) for the gap check
    tx2g = {}
    with open(os.path.join(toga, "orthology_classification.tsv")) as h:
        next(h)
        for line in h:
            f = line.rstrip("\n").split("\t")
            tx2g[f[1]] = f[0]
    gproj = collections.defaultdict(list)
    for name, loc in proj.items():
        tx = "#".join(name.split("#")[:2])
        if tx in tx2g:
            gproj[tx2g[tx]].append(loc)
    return status, o2o, gproj


def twobit_sizes(path):
    """Sequence sizes from a 2bit header (version 0 and 1; twoBitInfo is not installed everywhere)."""
    with open(path, "rb") as h:
        sig, ver, cnt, _ = struct.unpack("<IIII", h.read(16))
        if sig != 0x1A412743:
            sys.exit("not a little-endian 2bit file: " + path)
        off_fmt, off_len = ("<Q", 8) if ver == 1 else ("<I", 4)
        index = []
        for _ in range(cnt):
            n = h.read(1)[0]
            name = h.read(n).decode()
            index.append((name, struct.unpack(off_fmt, h.read(off_len))[0]))
        sizes = {}
        for name, off in index:
            h.seek(off)
            sizes[name] = struct.unpack("<I", h.read(4))[0]
    return sizes


def neighbours(gene, span, by_contig, status, o2o):
    c, s, e, st = span[gene]
    lst = by_contig[c]
    starts = [x[0] for x in lst]
    i = bisect.bisect_left(starts, s)
    left = right = None
    for j in range(i - 1, -1, -1):
        gs, ge, g = lst[j]
        if ge <= s and status.get(g) in INTACT and g in o2o:
            left = g
            break
    for j in range(i, len(lst)):
        gs, ge, g = lst[j]
        if gs >= e and status.get(g) in INTACT and g in o2o:
            right = g
            break
    return left, right


def wga_unaligned(paf, contigs):
    """Merged primary-aligned intervals of the reference (PAF query = reference genome,
    as in a reference-vs-query whole-genome minimap2 run) for the given contigs."""
    iv = collections.defaultdict(list)
    with open(paf) as h:
        for line in h:
            f = line.split("\t", 13)
            if f[0] in contigs and "tp:A:P" in line:
                iv[f[0]].append((int(f[2]), int(f[3])))
    merged = {}
    for c, xs in iv.items():
        xs.sort()
        m = []
        for s, e in xs:
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        merged[c] = m
    return merged


def unaligned_run(merged, c, s, e):
    """The reference stretch around the gene midpoint with no whole-genome alignment:
    (left, right), or None when the midpoint is aligned or the stretch reaches a contig end."""
    m = merged.get(c, [])
    mid = (s + e) // 2
    i = bisect.bisect_right([x[0] for x in m], mid) - 1
    if i >= 0 and m[i][1] > mid:
        return None
    if i < 0 or i + 1 >= len(m):
        return None
    return m[i][1], m[i + 1][0]


def gene_call(r):
    """One call per gene from all scales. An assembly gap anywhere between anchors blocks
    a deletion (EPS8L1/TMED4: a 100-N gap next to the site); a deletion needs an empty
    site and no scale that sees sequence in its place."""
    flank = [r.get(f"f{fl}_verdict", "") for fl in FLANKS] + [r.get("edge_verdict", "")]
    if "assembly_gap" in flank or r.get("nb_verdict") == "assembly_gap":
        return "assembly_gap"
    # only the fixed flanks can see the gene; edge "not_empty" is unaligned flank sequence
    # (toy: 6 of 100 cut genes became "conflict" when it counted)
    seen = {"present", "size_changed"}
    if "empty_site" in flank:
        return "conflict" if seen & set(flank[:len(FLANKS)]) else "deletion"
    if "present" in flank:
        return "present"
    return "no_call"


def pair_verdict(lhit, rhit, ref_gap, gene_len, qsizes, nrun, slack=0):
    """lhit/rhit: (contig, inner_pos, strand) or None. ref_gap = what lies
    between the two anchors in the reference (incl. the gene). slack = flank
    bases between the aligned pieces and the gene (extrapolated, not seen)."""
    if lhit is None and rhit is None:
        return "no_anchor", ""
    if lhit is None or rhit is None:
        c, p, st = (lhit or rhit)[:3]
        room = min(p, qsizes.get(c, 0) - p)
        return ("assembly_edge" if room < ref_gap else "one_flank"), ""
    if lhit[0] != rhit[0] or lhit[2] != rhit[2]:
        return "synteny_broken", ""
    gap = (rhit[1] - lhit[1]) if lhit[2] == "+" else (lhit[1] - rhit[1])
    # flanks may overlap in the query when the deletion took a little more than the gene
    if gap < -max(1000, EMPTY_FRAC * gene_len):
        return "synteny_broken", str(gap)
    # any N run between the anchors: an assembly gap of unknown size (scaffolders write 100 N),
    # so the interval length says nothing (b2t: 26 of 28 "large deletions" were one 100-N gap)
    # checked between the ALIGNED piece edges (element 3 when given), not the extrapolated
    # gene ends: EPS8L1/TMED4 were false empty sites, a flank extrapolated across a 100-N gap
    pts = [h[3] if len(h) > 3 else h[1] for h in (lhit, rhit)] + [lhit[1], rhit[1]]
    if nrun(lhit[0], min(pts), max(pts)) >= MIN_N_RUN:
        return "assembly_gap", str(gap)
    expected_empty = ref_gap - gene_len
    if gap - expected_empty <= EMPTY_FRAC * gene_len:
        if slack > max(MAX_SLACK, EMPTY_FRAC * gene_len):
            return "unresolved", str(gap)
        return "empty_site", str(gap)
    r = gap / ref_gap if ref_gap else 0
    if PRESENT_LO <= r <= PRESENT_HI:
        return "present", str(gap)
    return "size_changed", str(gap)


def main():
    ap = argparse.ArgumentParser()
    for a in ("toga", "ref-bed", "isoforms", "ref-2bit", "query-2bit", "query-index", "genes", "out"):
        ap.add_argument("--" + a, required=True)
    ap.add_argument("--threads", default="16")
    ap.add_argument("--no-neighbours", action="store_true", help="toy queries: skip the neighbour-gene scale")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--mm2-opts", default="-x asm20 -N 5")
    ap.add_argument("--wga-paf", help="reference-vs-query whole-genome PAF (reference = PAF query); adds the "
                    "'edge' scale: 1 kb flanks at the edges of the unaligned stretch around the gene, for "
                    "losses larger than the fixed flanks")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    span = read_genes(a.ref_bed, a.isoforms)
    status, o2o, gproj = read_toga(a.toga)
    by_contig = collections.defaultdict(list)
    for g, (c, s, e, st) in span.items():
        by_contig[c].append((s, e, g))
    for c in by_contig:
        by_contig[c].sort()
    rsizes = twobit_sizes(a.ref_2bit)
    qsizes = twobit_sizes(a.query_2bit)

    genes = [l.split()[0] for l in open(a.genes) if l.strip()]
    genes = [g for g in genes if g in span]

    # edge scale: unaligned stretch around the gene, only when it is longer than the
    # gene and shorter than MAX_EDGE_RUN (a whole unaligned chromosome arm is no site)
    runs = {}
    if a.wga_paf:
        merged = wga_unaligned(a.wga_paf, set(span[g][0] for g in genes))
        for g in genes:
            c, s, e, st = span[g]
            r = unaligned_run(merged, c, s, e)
            if r and r[1] - r[0] > e - s and r[1] - r[0] <= MAX_EDGE_RUN:
                runs[g] = r

    # flank BED: name = gene|side|size
    bed = os.path.join(a.out, "flanks.bed")
    with open(bed, "w") as h:
        for g in genes:
            c, s, e, st = span[g]
            if g in runs:
                rl, rr = runs[g]
                if rl >= EDGE_FLANK and rr + EDGE_FLANK <= rsizes[c]:
                    h.write(f"{c}\t{rl - EDGE_FLANK}\t{rl}\t{g}|L|edge\n{c}\t{rr}\t{rr + EDGE_FLANK}\t{g}|R|edge\n")
            for fl in FLANKS:
                ls, re_ = max(0, s - fl), min(rsizes[c], e + fl)
                if s - ls >= 100:
                    h.write(f"{c}\t{ls}\t{s}\t{g}|L|{fl}\n")
                if re_ - e >= 100:
                    h.write(f"{c}\t{e}\t{re_}\t{g}|R|{fl}\n")
    fa = os.path.join(a.out, "flanks.fa")
    paf = os.path.join(a.out, "flanks.paf")
    sh(f"twoBitToFa -bed={bed} {a.ref_2bit} {fa}")
    sh(f"{a.minimap2} -c {a.mm2_opts} -t {a.threads} {a.query_index} {fa} > {paf} 2> {paf}.log")

    # Pieces of one flank: minimap2 splits a flank at large indels, so keep all
    # primary pieces, choose the target (contig, strand) with the most aligned
    # flank bases, and place the gene end from the piece nearest the gene side.
    pieces = collections.defaultdict(list)
    with open(paf) as h:
        for line in h:
            f = line.split("\t")
            if "tp:A:P" not in line or int(f[11]) < MIN_MAPQ:
                continue
            tags = dict(t.split(":", 2)[::2] for t in f[12:])
            s1, s2 = int(tags.get("s1", 0)), int(tags.get("s2", 0))
            pieces[f[0]].append((int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8]),
                                 s2 / s1 if s1 else 0.0, s1))
    hits = {}
    for qn, ps in pieces.items():
        g, side, fl = qn.split("|")
        ql = ps[0][0]
        cov = collections.Counter()
        for _, qs, qe, strand, tn, ts, te, _, _ in ps:
            cov[(tn, strand)] += qe - qs
        (tn, strand), aligned = cov.most_common(1)[0]
        if aligned < min(MIN_QCOV * ql, MIN_ALIGNED):
            continue
        gpos = ql if side == "L" else 0          # gene side of the flank in flank coordinates
        best = min((p for p in ps if (p[4], p[3]) == (tn, strand)),
                   key=lambda p: abs((p[2] if side == "L" else p[1]) - gpos))
        _, qs, qe, _, _, ts, te, _, _ = best
        if strand == "+":
            inner = te + (gpos - qe) if side == "L" else ts - (qs - gpos)
        else:
            inner = ts - (gpos - qe) if side == "L" else te + (qs - gpos)
        edge = best[2] if side == "L" else best[1]
        # uniqueness from the strongest piece of the chosen target
        s2s1 = max((p for p in ps if (p[4], p[3]) == (tn, strand)), key=lambda p: p[8])[7]
        # target position of the aligned piece's gene-side end
        tedge = (te if side == "L" else ts) if strand == "+" else (ts if side == "L" else te)
        hits[(g, side, fl)] = (tn, inner, strand, abs(edge - gpos), s2s1, tedge)

    ncache = {}

    def nrun(c, s, e):
        """Longest N run in the query interval."""
        if e - s <= 0 or e - s > 5_000_000:
            return 0
        k = (c, s, e)
        if k not in ncache:
            seq = subprocess.run(["twoBitToFa", f"{a.query_2bit}:{c}:{s}-{e}", "stdout"],
                                 capture_output=True, text=True, check=True).stdout
            seq = "".join(seq.splitlines()[1:])
            ncache[k] = max((len(m) for m in re.findall("N+", seq.upper())), default=0)
        return ncache[k]

    out = os.path.join(a.out, "verdicts.tsv")
    cols = ["gene", "toga_status", "ref_locus", "gene_len"]
    for fl in FLANKS:
        cols += [f"f{fl}_verdict", f"f{fl}_query_gap", f"f{fl}_query"]
    cols += ["edge_run", "edge_verdict", "edge_query_gap", "edge_query"]
    cols += ["nb_left", "nb_right", "nb_verdict", "nb_query_gap", "nb_ref_gap", "query_projection_in_gap", "call"]
    with open(out, "w") as h:
        h.write("\t".join(cols) + "\n")
        for g in genes:
            c, s, e, st = span[g]
            glen = e - s
            row = [g, status.get(g, "NA"), f"{c}:{s}-{e}{st}", str(glen)]
            qint = None
            for fl in FLANKS:
                lh = hits.get((g, "L", str(fl)))
                rh = hits.get((g, "R", str(fl)))
                slack = (lh[3] if lh else 0) + (rh[3] if rh else 0)
                v, gap = pair_verdict(lh and (lh[0], lh[1], lh[2], lh[5]), rh and (rh[0], rh[1], rh[2], rh[5]),
                                      glen, glen, qsizes, nrun, slack)
                if any(x and x[4] >= NONUNIQUE_S2S1 for x in (lh, rh)):
                    v = "flank_not_unique"
                loc = ""
                if lh and rh and lh[0] == rh[0]:
                    lo, hi = sorted((lh[1], rh[1]))
                    loc = f"{lh[0]}:{lo}-{hi}{lh[2]}"
                    if qint is None:
                        qint = (lh[0], lo, hi)
                row += [v, gap, loc]
            ev = ["", "", "", ""]
            if g in runs:
                rl, rr = runs[g]
                lh, rh = hits.get((g, "L", "edge")), hits.get((g, "R", "edge"))
                slack = (lh[3] if lh else 0) + (rh[3] if rh else 0)
                v, gap = pair_verdict(lh and (lh[0], lh[1], lh[2], lh[5]), rh and (rh[0], rh[1], rh[2], rh[5]),
                                      rr - rl, rr - rl, qsizes, nrun, slack)
                # the stretch did not align, so query sequence in its place is not the gene:
                # only "empty or not" is meaningful at this scale (toy: unaligned flank parts
                # left 2.8-20 kb behind a gene-only cut)
                if v in ("present", "size_changed"):
                    v = "not_empty"
                if any(x and x[4] >= NONUNIQUE_S2S1 for x in (lh, rh)):
                    v = "flank_not_unique"
                loc = ""
                if lh and rh and lh[0] == rh[0]:
                    lo, hi = sorted((lh[1], rh[1]))
                    loc = f"{lh[0]}:{lo}-{hi}{lh[2]}"
                ev = [f"{c}:{rl}-{rr}", v, gap, loc]
            row += ev
            ln, rn = neighbours(g, span, by_contig, status, o2o)
            nv, ngap, rgap = "no_neighbour", "", ""
            if ln and rn and not a.no_neighbours:
                lq, rq = o2o[ln][1], o2o[rn][1]
                rgap_i = span[rn][1] - span[ln][2]
                same = span[ln][3] == lq[3]          # neighbour keeps its strand -> '+' orientation
                lhit = (lq[0], lq[2] if same else lq[1], "+" if same else "-")
                rhit = (rq[0], rq[1] if same else rq[2], "+" if same else "-")
                # neighbours are tens of kb apart: intergenic size change hides
                # a gene-sized deletion, so this scale tests synteny only
                nv, ngap = pair_verdict(lhit, rhit, rgap_i, 0, qsizes, nrun)
                if nv in ("empty_site", "present", "size_changed"):
                    r = int(ngap) / rgap_i if rgap_i else 0
                    nv = "syntenic" if PRESENT_LO <= r <= PRESENT_HI else "syntenic_size_changed"
                rgap = str(rgap_i)
            row += [ln or "", rn or "", nv, ngap, rgap]
            inside = ""
            if qint:
                inside = ",".join(f"{p[0]}:{p[1]}-{p[2]}" for p in gproj.get(g, [])
                                  if p[0] == qint[0] and p[1] < qint[2] and p[2] > qint[1])
            row.append(inside)
            row.append(gene_call(dict(zip(cols, row))))
            h.write("\t".join(row) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
