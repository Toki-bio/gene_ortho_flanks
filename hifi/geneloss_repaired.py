#!/usr/bin/env python3
"""Re-call betulina -> trizona gene losses on BOTH trizona haplotypes (gap-repaired deposited hap1/hap2, and their
gene-search versions), with gof.py's flank-anchor logic plus an exon check.

Per gene (reference = betulina, gof.py b2t inputs): the fixed flanks (1 kb, 10 kb) and the edge flanks of the
unaligned stretch (whole-genome PAF, gof.py --wga-paf) are placed in one query assembly with gof.place_flanks, each
scale is judged with gof.pair_verdict (N run >= 10 between aligned anchor edges -> assembly_gap); the gene region
(gene +- 10 kb) is aligned too (--cs) and the union of the gene's exons (all its transcripts in the BED12) is
counted as aligned bases on the anchor locus.

Per haplotype call:
  present       exons >= 80 % aligned at the anchored locus (or, with no anchor, at the best primary locus:
                present_unanchored)
  present_diverged  20-80 % of exon bases aligned, a fixed-flank scale says present (locus interval 0.5-2x the
                gene) and no N run in the window: the locus is intact, exons partly too diverged for asm20
  deleted       a scale says empty_site, no N run in the locus window, exons < 20 % aligned
  partial       20-80 % of exon bases aligned, no N run in the window, no fixed scale present (part missing?)
  gap           an N run in the locus window (or assembly_gap at any scale) and exons < 80 % aligned; the gap is
                identified by the 40 bp left of its N run (original deposited gaps.tsv) -> gap id + cause
  absent_chrom  X-linked gene (anchors on dep1 OZ418353) on hap2, which has no X
  no_call       none of the above (no anchors, flanks not unique, synteny broken ...)

Subcommands:
  prep   VERDICTS REF_BED ISOFORMS REF_2BIT WGA_PAF OUT [EXTRA_GENES]
         genes with call deletion/conflict/assembly_gap (+ EXTRA_GENES: comma list or @file, e.g. FI controls) -> OUT/genes.tsv, OUT/query.bed,
         OUT/query.fa (twoBitToFa), OUT/exons.tsv
  call   OUT LABEL PAF ASM.fa ORIG.fa GAPS.tsv CAUSES.tsv  -> OUT/LABEL.calls.tsv
  combine OUT  (needs dep1, dep2, dep1gs, dep2gs calls) -> OUT/final.tsv, OUT/summary.txt
"""
import collections
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import gof  # noqa: E402

REGION_FLANK = 10000
PRESENT_FRAC, DELETED_FRAC = 0.8, 0.2
WIN_PAD = 1000            # locus window = anchors +- this
LIST_PAD = 20000          # N runs this close to the window are listed (not called) as nearby gaps
GAP_KEY = 40              # bases left of an N run that identify the deposited gap
TARGET_CALLS = ("deletion", "conflict", "assembly_gap")


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def faidx(fa, regions):
    """{region: sequence} via samtools faidx -r (regions 'chrom:start-end', 1-based inclusive)."""
    if not regions:
        return {}
    rf = f"{os.environ.get('TMPDIR', '/tmp')}/gl2_reg_{os.getpid()}.txt"
    with open(rf, "w") as h:
        h.write("\n".join(regions) + "\n")
    out = subprocess.run(["samtools", "faidx", "-r", rf, fa], capture_output=True, text=True, check=True).stdout
    os.remove(rf)
    seqs, name = {}, None
    for l in out.splitlines():
        if l.startswith(">"):
            name = l[1:]
            seqs[name] = []
        else:
            seqs[name].append(l)
    return {k: "".join(v) for k, v in seqs.items()}


def read_fai(fa):
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fa + ".fai")}


# ---------------------------------------------------------------- prep
def prep(ver, ref_bed, isoforms, ref2bit, wga_paf, out, extra=""):
    os.makedirs(out, exist_ok=True)
    rows = [l.rstrip("\n").split("\t") for l in open(ver)]
    hd, rows = rows[0], rows[1:]
    ix = {k: i for i, k in enumerate(hd)}
    if extra.startswith("@"):
        extra = ",".join(l.split()[0] for l in open(extra[1:]) if l.strip())
    extra = set(x for x in extra.split(",") if x)
    sel = [r for r in rows if r[ix["call"]] in TARGET_CALLS or r[0] in extra]
    span = gof.read_genes(ref_bed, isoforms)
    rsizes = gof.twobit_sizes(ref2bit)
    genes = [r[0] for r in sel if r[0] in span]
    # exons: union over all transcripts of the gene
    tx2gene = {l.split("\t")[1].strip(): l.split("\t")[0] for l in open(isoforms)}
    ex = collections.defaultdict(list)
    for l in open(ref_bed):
        f = l.rstrip("\n").split("\t")
        g = tx2gene.get(f[3])
        if g not in span or f[0] != span[g][0]:
            continue
        st = int(f[1])
        for sz, of in zip(f[10].strip(",").split(","), f[11].strip(",").split(",")):
            ex[g].append((st + int(of), st + int(of) + int(sz)))
    merged = gof.wga_unaligned(wga_paf, set(span[g][0] for g in genes))
    with open(f"{out}/genes.tsv", "w") as hg, open(f"{out}/query.bed", "w") as hb, open(f"{out}/exons.tsv", "w") as he:
        hg.write("gene\ttoga_status\tprev_call\tref_locus\tgene_len\tregion\n")
        for r in sel:
            g = r[0]
            if g not in span:
                print("not in BED:", g, file=sys.stderr)
                continue
            c, s, e, st = span[g]
            rs, re_ = max(0, s - REGION_FLANK), min(rsizes[c], e + REGION_FLANK)
            hg.write(f"{g}\t{r[ix['toga_status']]}\t{r[ix['call']]}\t{c}:{s}-{e}{st}\t{e - s}\t{c}:{rs}-{re_}\n")
            hb.write(f"{c}\t{rs}\t{re_}\t{g}|reg|{REGION_FLANK}\n")
            run = gof.unaligned_run(merged, c, s, e)
            if run and run[1] - run[0] > e - s and run[1] - run[0] <= gof.MAX_EDGE_RUN:
                rl, rr = run
                if rl >= gof.EDGE_FLANK and rr + gof.EDGE_FLANK <= rsizes[c]:
                    hb.write(f"{c}\t{rl - gof.EDGE_FLANK}\t{rl}\t{g}|L|edge\n{c}\t{rr}\t{rr + gof.EDGE_FLANK}\t{g}|R|edge\n")
            for fl in gof.FLANKS:
                ls, rr2 = max(0, s - fl), min(rsizes[c], e + fl)
                if s - ls >= 100:
                    hb.write(f"{c}\t{ls}\t{s}\t{g}|L|{fl}\n")
                if rr2 - e >= 100:
                    hb.write(f"{c}\t{e}\t{rr2}\t{g}|R|{fl}\n")
            xs = sorted(set(ex[g]))
            m = []
            for a, b in xs:
                if m and a <= m[-1][1]:
                    m[-1][1] = max(m[-1][1], b)
                else:
                    m.append([a, b])
            # exon coordinates relative to the region query (forward reference strand)
            he.write(g + "\t" + ",".join(f"{a - rs}-{b - rs}" for a, b in m) + "\n")
    sh(f"twoBitToFa -bed={out}/query.bed {ref2bit} {out}/query.fa")
    print(f"prep: {len(genes)} genes", file=sys.stderr)


# ---------------------------------------------------------------- call
def cs_aligned(qs, strand, ql, cs):
    """Query intervals (forward query coordinates) aligned (match or mismatch) in one PAF record."""
    out, q = [], 0          # q = offset along the alignment in query-strand-of-alignment order
    for op, val in re.findall(r"([:*+\-~])([0-9]+|[a-z]+|[a-z]{2}[0-9]+[a-z]{2})", cs):
        if op == ":":
            n = int(val)
            out.append((q, q + n))
            q += n
        elif op == "*":
            out.append((q, q + 1))
            q += 1
        elif op == "+":
            q += len(val)
        # '-' (target-only bases) and '~' (introns) do not consume query
    span = []
    alen = q
    for a, b in out:
        if strand == "+":
            span.append((qs + a, qs + b))
        else:
            qe = qs + alen
            span.append((qe - b, qe - a))
    return span


def covered(exons, ivs):
    """(exon bp, aligned exon bp, exons >= 50 % aligned); overlapping alignment pieces are merged first."""
    m = []
    for x, y in sorted(ivs):
        if m and x <= m[-1][1]:
            m[-1][1] = max(m[-1][1], y)
        else:
            m.append([x, y])
    tot = cov = nex = 0
    for a, b in exons:
        c = sum(max(0, min(b, y) - max(a, x)) for x, y in m)
        tot += b - a
        cov += c
        nex += c >= 0.5 * (b - a)
    return tot, cov, nex


def call(out, label, paf, asm, orig, gaps_tsv, causes_tsv):
    genes = [l.rstrip("\n").split("\t") for l in open(f"{out}/genes.tsv")][1:]
    exons = {}
    for l in open(f"{out}/exons.tsv"):
        g, s = l.rstrip("\n").split("\t")
        exons[g] = [tuple(map(int, x.split("-"))) for x in s.split(",") if x]
    sizes = read_fai(asm)
    hits = gof.place_flanks(paf)
    # region alignments (all records, primary and secondary)
    reg = collections.defaultdict(list)
    for l in open(paf):
        f = l.rstrip("\n").split("\t")
        g, side, _ = f[0].split("|")
        if side != "reg":
            continue
        cs = next((t[5:] for t in f[12:] if t.startswith("cs:Z:")), "")
        tp = next((t[5:] for t in f[12:] if t.startswith("tp:A:")), "")
        reg[g].append(dict(ql=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], t=f[5], ts=int(f[7]),
                           te=int(f[8]), nm=int(f[9]), bl=int(f[10]), mapq=int(f[11]), tp=tp,
                           ivs=cs_aligned(int(f[2]), f[4], int(f[1]), cs)))
    # deposited gaps -> key (GAP_KEY bases left of the N run) -> gap id, cause
    cause = {l.split("\t")[0]: l.rstrip("\n").split("\t")[-1] for l in list(open(causes_tsv))[1:]}
    gl = [l.rstrip("\n").split("\t") for l in list(open(gaps_tsv))[1:]]
    regs = {f"{f[1]}:{int(f[2]) - GAP_KEY + 1}-{int(f[2])}": f[0] for f in gl if int(f[2]) > GAP_KEY}
    keyseq = faidx(orig, list(regs))
    gapkey = {}
    for r, gid in regs.items():
        gapkey.setdefault(keyseq[r].upper(), []).append(gid)

    seq_cache = {}

    def window_seq(c, s, e):
        s, e = max(0, s), min(sizes.get(c, 0), e)
        if e <= s:
            return s, ""
        k = (c, s, e)
        if k not in seq_cache:
            seq_cache[k] = faidx(asm, [f"{c}:{s + 1}-{e}"]).get(f"{c}:{s + 1}-{e}", "").upper()
        return s, seq_cache[k]

    def nruns(c, s, e):
        """N runs >= MIN_N_RUN in [s, e): list of (start, end, gap_id|?, cause)."""
        if e - s > 5_000_000:
            return []
        off, sq = window_seq(c, s - GAP_KEY, e)
        res = []
        for m in re.finditer("N{%d,}" % gof.MIN_N_RUN, sq):
            a, b = off + m.start(), off + m.end()
            if b <= s:
                continue
            left = sq[max(0, m.start() - GAP_KEY):m.start()]
            ids = gapkey.get(left, []) if len(left) == GAP_KEY else []
            gid = ids[0] if len(ids) == 1 else ("?" if not ids else "|".join(ids))
            res.append((a, b, gid, cause.get(gid, "unknown")))
        return res

    def nrun_len(c, s, e):
        return max((b - a for a, b, _, _ in nruns(c, s, e)), default=0)

    cols = ["gene", "toga_status", "prev_call", "ref_locus", "gene_len", "scale_verdicts", "anchor_scale", "locus",
            "exon_bp", "exon_aligned", "exons_covered", "exon_frac", "exon_ident", "n_runs_in_window", "gaps_in_window",
            "gaps_nearby", "call"]
    rows = []
    for g, tstat, prev, rloc, glen, region in genes:
        glen = int(glen)
        verd, pairs = {}, {}
        for fl in [str(x) for x in gof.FLANKS]:
            lh, rh = hits.get((g, "L", fl)), hits.get((g, "R", fl))
            slack = (lh[3] if lh else 0) + (rh[3] if rh else 0)
            v, _ = gof.pair_verdict(lh and (lh[0], lh[1], lh[2], lh[5]), rh and (rh[0], rh[1], rh[2], rh[5]),
                                    glen, glen, sizes, nrun_len, slack)
            if any(x and x[4] >= gof.NONUNIQUE_S2S1 for x in (lh, rh)):
                v = "flank_not_unique"
            verd[fl] = v
            pairs[fl] = (lh, rh)
        pairs["edge"] = (hits.get((g, "L", "edge")), hits.get((g, "R", "edge")))
        # edge scale: needs the reference stretch length = distance of the edge flanks in the query.bed
        lh, rh = pairs["edge"]
        if lh or rh:
            ev, _ = gof.pair_verdict(lh and (lh[0], lh[1], lh[2], lh[5]), rh and (rh[0], rh[1], rh[2], rh[5]),
                                     EDGE_LEN[g], EDGE_LEN[g], sizes, nrun_len,
                                     (lh[3] if lh else 0) + (rh[3] if rh else 0))
            if ev in ("present", "size_changed"):
                ev = "not_empty"
            if any(x and x[4] >= gof.NONUNIQUE_S2S1 for x in (lh, rh)):
                ev = "flank_not_unique"
            verd["edge"] = ev
        else:
            verd["edge"] = "no_anchor" if g in EDGE_LEN else ""
        # locus window: first scale with both anchors on one contig/strand, else one anchor
        win, ascale = None, ""
        for fl in ["10000", "1000", "edge"]:
            lh, rh = pairs[fl]
            if lh and rh and lh[0] == rh[0] and lh[2] == rh[2] and verd[fl] != "flank_not_unique":
                pts = [lh[1], rh[1], lh[5], rh[5]]
                if max(pts) - min(pts) > 3 * (EDGE_LEN.get(g, 0) + glen + 2 * REGION_FLANK) + 100_000:
                    continue          # anchors far apart: not one locus (synteny broken)
                win = (lh[0], min(pts) - WIN_PAD, max(pts) + WIN_PAD, lh[2])
                ascale = fl + ":pair"
                break
        if win is None:
            for fl in ["10000", "1000", "edge"]:
                lh, rh = pairs[fl]
                h = lh or rh
                if h and h[4] < gof.NONUNIQUE_S2S1:
                    ext = (EDGE_LEN.get(g, 0) if fl == "edge" else 0) + glen + 2 * REGION_FLANK
                    # extend towards the gene side
                    toward_plus = (h is lh) == (h[2] == "+")
                    win = (h[0], h[1] - WIN_PAD, h[1] + ext) if toward_plus else (h[0], h[1] - ext, h[1] + WIN_PAD)
                    win = win + (h[2],)
                    ascale = fl + ":one"
                    break
        # exon coverage on the window (or the best primary locus when unanchored)
        recs = reg.get(g, [])
        if win:
            use = [r for r in recs if r["t"] == win[0] and r["strand"] == win[3] and r["te"] > win[1] and r["ts"] < win[2]]
        else:
            prim = [r for r in recs if r["tp"] == "P"]
            best = collections.Counter()
            for r in prim:
                best[(r["t"], r["strand"])] += r["qe"] - r["qs"]
            use = []
            if best:
                (t, st), _ = best.most_common(1)[0]
                use = [r for r in prim if (r["t"], r["strand"]) == (t, st)]
                if use:
                    ascale = "unanchored"
                    win = (t, min(r["ts"] for r in use) - WIN_PAD, max(r["te"] for r in use) + WIN_PAD, st)
        ivs = [iv for r in use for iv in r["ivs"]]
        tot, cov, nex = covered(exons.get(g, []), ivs)
        frac = cov / tot if tot else 0.0
        bl = sum(r["bl"] for r in use)
        ident = (sum(r["nm"] for r in use) / bl) if bl else 0.0   # PAF col 10 = matching bases
        inw, near = [], []
        if win:
            c, s, e = win[0], max(0, win[1]), min(sizes.get(win[0], 0), win[2])
            inw = nruns(c, s, e)
            near = [x for x in nruns(c, s - LIST_PAD, e + LIST_PAD) if x not in inw]
        sv = set(verd.values())
        if frac >= PRESENT_FRAC:
            res = "present" if ascale != "unanchored" else "present_unanchored"
        elif inw or "assembly_gap" in sv:
            res = "gap"
        elif "empty_site" in sv and frac < DELETED_FRAC:
            res = "deleted"
        elif frac >= DELETED_FRAC and win and ascale != "unanchored" and                 any(verd[k] == "present" for k in ("1000", "10000")):
            # locus intact at a fixed-flank scale (interval 0.5-2x the gene), no N: exons only partly aligned
            # because of divergence (FI control USP22: 10/13 exons, 0.755 on hap2)
            res = "present_diverged"
        elif frac >= DELETED_FRAC and win and ascale != "unanchored":
            res = "partial"
        else:
            res = "no_call"
        locus = f"{win[0]}:{max(0, win[1])}-{win[2]}{win[3]}" if win else ""
        fmt = lambda xs: ",".join(f"{gid}:{cz}" for _, _, gid, cz in xs)
        rows.append([g, tstat, prev, rloc, str(glen),
                     ";".join(f"{k}={verd[k] or '-'}" for k in ["1000", "10000", "edge"]), ascale, locus,
                     str(tot), str(cov), f"{nex}/{len(exons.get(g, []))}", f"{frac:.3f}", f"{ident:.4f}", str(len(inw)), fmt(inw), fmt(near), res])
    with open(f"{out}/{label}.calls.tsv", "w") as h:
        h.write("\t".join(cols) + "\n")
        for r in rows:
            h.write("\t".join(r) + "\n")
    print(label, collections.Counter(r[-1] for r in rows), file=sys.stderr)


EDGE_LEN = {}


def load_edge_len(out):
    """reference stretch between edge flanks (query.bed) per gene, for the edge-scale pair_verdict."""
    e = {}
    for l in open(f"{out}/query.bed"):
        c, s, en, n = l.rstrip("\n").split("\t")
        g, side, fl = n.split("|")
        if fl == "edge":
            e.setdefault(g, {})[side] = (int(s), int(en))
    for g, d in e.items():
        if "L" in d and "R" in d:
            EDGE_LEN[g] = d["R"][0] - d["L"][1]


# ---------------------------------------------------------------- combine
def combine(out, xchrom="OZ418353"):
    def load(lb):
        p = f"{out}/{lb}.calls.tsv"
        if not os.path.exists(p):
            return {}
        rr = [l.rstrip("\n").split("\t") for l in open(p)]
        return {r[0]: dict(zip(rr[0], r)) for r in rr[1:]}
    d1, d2, g1, g2 = load("dep1"), load("dep2"), load("dep1gs"), load("dep2gs")
    o1, o2 = load("orig1"), load("orig2")      # unrepaired deposited haplotypes (before repair)
    gen = [l.rstrip("\n").split("\t") for l in open(f"{out}/genes.tsv")][1:]
    tally, lists = collections.Counter(), collections.defaultdict(list)
    prevtally = collections.Counter()
    gapcause = collections.Counter()

    def short(r):
        if r["call"] == "gap":
            return "gap(" + (r["gaps_in_window"] or "assembly_gap_scale") + ")"
        return r["call"]

    with open(f"{out}/final.tsv", "w") as h:
        h.write("gene\ttoga_status\tprev_call\tref_locus\tgene_len\thap1\thap2\thap1_exon_frac\thap2_exon_frac\t"
                "hap1_exons\thap2_exons\thap1_locus\thap2_locus\tanimal\thap1_genesearch\thap2_genesearch\t"
                "readable_copy_note\thap1_before_repair\thap2_before_repair\n")
        for g, ts, prev, rloc, glen, _ in gen:
            a, b = d1[g], d2[g]
            ca, cb = a["call"], b["call"]
            xl = a["locus"].startswith(xchrom)
            P = ("present", "present_unanchored", "present_diverged")
            if xl and cb not in P:
                cb = "absent_chrom"       # hap2 has no X (male): hap1 is the only copy
            if xl:
                animal = {"present": "present (X, hemizygous)", "present_unanchored": "present (X, hemizygous, unanchored)", "present_diverged": "present (X, hemizygous, diverged)",
                          "deleted": "deletion (X, hemizygous)", "partial": "partial deletion (X, hemizygous)",
                          "gap": "unknown (X, gap)"}.get(ca, "unknown (X, no_call)")
            elif ca in P and cb in P:
                animal = "present_both"
            elif {ca, cb} == {"deleted"}:
                animal = "homozygous_deletion"
            elif (ca in P and cb == "deleted") or (cb in P and ca == "deleted"):
                animal = "heterozygous_deletion"
            elif (ca in P and cb == "partial") or (cb in P and ca == "partial"):
                animal = "heterozygous_partial_deletion"
            elif ca == "partial" and cb == "partial":
                animal = "homozygous_partial_deletion"
            elif {ca, cb} <= {"deleted", "partial"}:
                animal = "homozygous_deletion_or_partial"
            elif ca in P or cb in P:
                other = cb if ca in P else ca
                animal = f"present_one_hap (other: {other})"
            elif "gap" in (ca, cb) and {ca, cb} & {"deleted", "partial"}:
                animal = "unknown_gap_other_deleted"
            elif "gap" in (ca, cb):
                animal = "unknown_gap"
            else:
                animal = "unknown_no_call"
            note = ""
            ga, gb = g1.get(g, {}).get("call", ""), g2.get(g, {}).get("call", "")
            PG = ("present", "present_unanchored", "present_diverged")
            if ca not in P and ga in PG:
                note += "hap1 present on readable copy only;"
            if cb not in P and cb != "absent_chrom" and gb in PG:
                note += "hap2 present on readable copy only;"
            h.write("\t".join([g, ts, prev, rloc, glen, short(a),
                               cb if cb == "absent_chrom" else short(b), a["exon_frac"], b["exon_frac"],
                               a["exons_covered"], b["exons_covered"], a["locus"], b["locus"], animal, ga, gb, note,
                               o1.get(g, {}).get("call", ""), o2.get(g, {}).get("call", "")]) + "\n")
            if "gap" in (ca, cb):
                for r in (a, b):
                    if r["call"] != "gap":
                        continue
                    xs = [x for x in r["gaps_in_window"].split(",") if x] or ["-:N outside window (assembly_gap scale)"]
                    for x in xs:
                        gapcause[x.split(":", 1)[1]] += 1
            key = animal.split(" (")[0]
            tally[key] += 1
            prevtally[(prev, key)] += 1
            lists[key].append(g)
    with open(f"{out}/summary.txt", "w") as h:
        h.write("animal-level call\tn\n")
        for k, n in tally.most_common():
            h.write(f"{k}\t{n}\n")
        h.write("\ngap causes of N runs in the locus window, haplotypes called gap (one count per haplotype x gap)\n")
        for k, n in gapcause.most_common():
            h.write(f"{k}\t{n}\n")
        h.write("\nprev_call x animal\n")
        for (p, k), n in sorted(prevtally.items()):
            h.write(f"{p}\t{k}\t{n}\n")
        for k in sorted(lists):
            if "deletion" in k or "unknown" in k or "present_one_hap" in k:
                h.write(f"\n{k}: {' '.join(sorted(lists[k]))}\n")
    print(open(f"{out}/summary.txt").read())


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "prep":
        prep(*sys.argv[2:])
    elif cmd == "call":
        load_edge_len(sys.argv[2])
        call(*sys.argv[2:])
    elif cmd == "combine":
        combine(*sys.argv[2:])
    else:
        sys.exit(__doc__)
