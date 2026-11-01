#!/usr/bin/env python3
"""Copy-number-free test of whether a repaired gap piece is carried by BOTH chromosome copies (category 1, hom_one_hap
inserts) or by one only (category 2, one-copy-unreadable genesearch inserts). No read depth, no k-mer counts are used
for the call; only which haplotype's molecules continue into the piece.

Test A (HiFi, phased reads): heterozygous SNPs are called from the reads themselves in the flanks (+-F around each
junction, genome-wide mapping to the haploid gene-search assembly, so MAPQ is meaningful). A read alignment is phased
'ref' or 'alt' by the alleles it carries at flank het sites (>= 1 informative site, >= 80 % one allele). Per junction we
count phased reads that CROSS (one alignment reaching >= 1 kb into both flank and piece) and phased reads that come
from the flank and STOP within 1.5 kb of the junction (clip = >= 1 kb of read left unaligned; readend otherwise).
  both  : ref and alt each >= 2 crossing reads  -> two copies carry the piece
  one   : one allele >= 4 crossing, the other 0 -> one copy
Het sites INSIDE the piece (alt fraction 0.25-0.75) are counted too: they mean two different piece copies exist.

Test B (Hi-C, allele-specific contacts): Hi-C pairs fished by 31-mers carrying either allele of the flank het sites,
both mates mapped single-end (minimap2 -x sr) to the same haploid assembly. For each mate that covers a het site
(allele from cs) whose partner maps uniquely (MAPQ >= 20) in the same region, the partner is classed piece / same-side
flank (>= 1 kb from the site) / other-side flank. Two-copy prediction: the fraction piece/(piece+same flank) is the same
for both alleles (identical distance distribution, because both alleles sit at the same sites). One-copy prediction:
the non-carrier allele reaches the piece only by trans-homolog contacts (near 0).

Subcommands
  prep     GS1.fa BED1 GS2.fa BED2 OUT [--ncat2 N --nctl N --hetctl NAME:FA:START:END ...]
           -> OUT/loci.tsv, OUT/regions.bed (asm chrom s e), OUT/regions.fa, OUT/hifi_cand.fa (25-mers every 40 bp)
  anchors  CAND.fa COUNTS(jellyfish dump -c) MINC MAXC > anchors.txt
  filter   REGIONS.bed ASM < PAF > filtered.paf    (primary, MAPQ >= 10, block >= 1 kb, identity >= 0.97, in regions)
  het      LOCI REGIONS FILT1 FILT2 > het.tsv
  phase    LOCI HET FILT1 FILT2 OUTPREFIX -> OUTPREFIX.junctions.tsv, OUTPREFIX.loci.tsv, OUTPREFIX.summary.txt
  hicanc   LOCI HET REGIONS.fa OUT -> OUT/hic_sites.tsv, OUT/hic_cand.fa (31-mers, both alleles)
  hic      LOCI HIC_SITES PAF1 PAF2 OUTPREFIX (PAF: minimap2 -x sr --cs of mates named READ_1 / READ_2)
  betq     LOCI GS1 GS2 MAXLEN > query.fa;  bet LOCI PAF(asm20 vs betulina) > bet.tsv"""
import bisect
import gzip
import math
import random
import re
import sys
from collections import Counter, defaultdict

F = 200000            # region half-width around each junction
J_IN = 1000           # crossing: alignment covers >= J_IN on both sides of the junction
STOP_WIN = 1500       # non-crossing flank alignment ending within this of the junction
APPROACH = 20000      # non-crossing flank alignments ending within this are counted as approaching (closest end reported)
NEAR = 10000          # flank het sites counted as 'near' the junction
MIN_DEPTH, MIN_ALT, LO, HI = 8, 3, 0.25, 0.75
K_HIC = 31
tok = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+)")
comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def opn(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_fa(path, upper=True):
    seqs, name, buf = {}, None, []
    with opn(path) as fh:
        for l in fh:
            if l.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = l[1:].split()[0], []
            else:
                buf.append(l.strip().upper() if upper else l.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


# ---------------------------------------------------------------- prep
def prep(a):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("gs1"); ap.add_argument("bed1"); ap.add_argument("gs2"); ap.add_argument("bed2"); ap.add_argument("out")
    ap.add_argument("--ncat2", type=int, default=60)
    ap.add_argument("--nctl", type=int, default=40)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--xchrom", default="OZ418353.1")
    ap.add_argument("--hetctl", action="append", default=[], help="NAME:FA:START:END piece inside FA (asm dep2)")
    ap.add_argument("--only", default="", help="comma list of locus ids (toy)")
    o = ap.parse_args(a)
    rnd = random.Random(o.seed)
    gs = {"dep1": read_fa(o.gs1), "dep2": read_fa(o.gs2)}
    loci = []
    beds = {}
    for asm, bed in (("dep1", o.bed1), ("dep2", o.bed2)):
        rows = [l.rstrip("\n").split("\t") for l in open(bed)]
        beds[asm] = rows
        cat2 = []
        for c, s, e, nm in rows:
            gap, cause, act = nm.split("|")
            s, e = int(s), int(e)
            if act == "insert":
                loci.append([f"{asm}_{gap}", "cat1", asm, c, s, e, cause])
            elif act == "genesearch" and e - s >= 2000:
                cat2.append([f"{asm}_{gap}", "cat2", asm, c, s, e, cause])
        rnd.shuffle(cat2)
        loci += sorted(cat2[:o.ncat2 // 2], key=lambda r: r[0])
    lens = [r[5] - r[4] for r in loci if r[1] == "cat1"]
    # random two-copy controls on dep1 autosomes, far from any N run and any repair
    g1 = gs["dep1"]
    auto = [c for c in g1 if c.startswith("OZ4183") and c != o.xchrom and len(g1[c]) > 20_000_000]
    bad = defaultdict(list)
    for c in auto:
        for m in re.finditer(r"N{10,}", g1[c]):
            bad[c].append((m.start(), m.end()))
    for c, s, e, nm in beds["dep1"]:
        bad[c].append((int(s), int(e)))
    for c in bad:
        bad[c].sort()
    tot = sum(len(g1[c]) for c in auto)
    n = 0
    while n < o.nctl:
        x = rnd.randrange(tot)
        for c in auto:
            if x < len(g1[c]):
                break
            x -= len(g1[c])
        L = min(100000, max(2000, rnd.choice(lens)))
        s, e = x, x + L
        lo, hi = s - F - 50000, e + F + 50000
        if lo < 0 or hi > len(g1[c]):
            continue
        if any(bs < hi and be > lo for bs, be in bad[c]):
            continue
        n += 1
        loci.append([f"ctl{n:02d}", "ctl_rand", "dep1", c, s, e, "random"])
    # known heterozygous deletions (one-copy controls) located in dep2 gene-search coordinates
    for h in o.hetctl:
        name, fa, ps, pe = h.split(":")
        w = read_fa(fa)
        (hdr, seq), = w.items()
        c = hdr.split(":")[0]
        pos = gs["dep2"][c].find(seq[:3000])
        if pos < 0 or gs["dep2"][c][pos:pos + len(seq)] != seq:
            print(f"WARN {name}: window not found exactly in dep2 {c}", file=sys.stderr)
            continue
        loci.append([name, "ctl_het", "dep2", c, pos + int(ps), pos + int(pe), "het_deletion"])
    if o.only:
        keep = set(o.only.split(","))
        loci = [r for r in loci if r[0] in keep]
    with open(f"{o.out}/loci.tsv", "w") as fh:
        fh.write("locus\tcat\tasm\tchrom\tps\tpe\tcause\n")
        for r in loci:
            fh.write("\t".join(map(str, r)) + "\n")
    reg = defaultdict(list)
    for lid, cat, asm, c, s, e, cause in loci:
        L = len(gs[asm][c])
        for j in (s, e):
            reg[(asm, c)].append((max(0, j - F), min(L, j + F)))
    with open(f"{o.out}/regions.bed", "w") as fb, open(f"{o.out}/regions.fa", "w") as ff, \
            open(f"{o.out}/hifi_cand.fa", "w") as fc:
        nk = 0
        for (asm, c), iv in sorted(reg.items()):
            for s, e in merge(iv):
                fb.write(f"{asm}\t{c}\t{s}\t{e}\n")
                seq = gs[asm][c][s:e]
                ff.write(f">{asm}|{c}:{s}-{e}\n{seq}\n")
                for i in range(0, len(seq) - 25, 40):
                    k = seq[i:i + 25]
                    if "N" not in k:
                        nk += 1
                        fc.write(f">{nk}\n{k}\n")
    print(Counter(r[1] for r in loci), file=sys.stderr)


def anchors(a):
    """k-mers of CAND whose count in COUNTS (k-mers missing from the dump count 0) is within [MINC, MAXC]"""
    cand, counts, minc, maxc = a[0], a[1], int(a[2]), int(a[3])
    cnt = {}
    for l in open(counts):
        k, n = l.split()
        cnt[k] = int(n)
    rc = lambda s: s.translate(comp)[::-1]
    seen = set()
    for l in open(cand):
        if l.startswith(">"):
            continue
        k = l.strip()
        cn = min(k, rc(k))
        n = cnt.get(k, cnt.get(rc(k), 0))
        if minc <= n <= maxc and cn not in seen:
            seen.add(cn)
            print(k)


# ---------------------------------------------------------------- filter
def load_regions(path):
    r = defaultdict(list)
    for l in open(path):
        asm, c, s, e = l.split()
        r[(asm, c)].append((int(s), int(e)))
    return r


def filt(a):
    reg = load_regions(a[0])
    asm = a[1]
    idx = {c: sorted(v) for (x, c), v in reg.items() if x == asm}
    starts = {c: [s for s, e in v] for c, v in idx.items()}
    for l in sys.stdin:
        f = l.rstrip("\n").split("\t")
        t = f[5]
        if t not in idx:
            continue
        tags = f[12:]
        if "tp:A:P" not in tags:
            continue
        ts, te, nm, bl, mq = int(f[7]), int(f[8]), int(f[9]), int(f[10]), int(f[11])
        if mq < 10 or bl < 1000 or nm / bl < 0.97:
            continue
        i = bisect.bisect_right(starts[t], te) - 1      # regions are merged: only region i can overlap
        if i < 0 or idx[t][i][1] <= ts:
            continue
        cs = next(x for x in tags if x.startswith("cs:Z:"))
        print("\t".join(f[:12] + [cs]))


def load_paf(path, asm):
    al = []
    with opn(path) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            al.append((f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8]), f[12][5:], asm))
    return al


def walk(ts, cs, sites):
    """alleles of one alignment at sorted ref positions `sites` (list) -> {pos: base}; '=' ref, '-' deleted."""
    out = {}
    if not sites:
        return out
    p = ts
    j = bisect.bisect_left(sites, ts)
    for m in tok.findall(cs):
        if j >= len(sites):
            break
        if m[0] == ":":
            n = int(m[1:])
            while j < len(sites) and sites[j] < p + n:
                if sites[j] >= p:
                    out[sites[j]] = "="
                j += 1
            p += n
        elif m[0] == "*":
            if sites[j] == p:
                out[p] = m[2].upper()
                j += 1
            p += 1
        elif m[0] == "-":
            n = len(m) - 1
            while j < len(sites) and sites[j] < p + n:
                if sites[j] >= p:
                    out[sites[j]] = "-"
                j += 1
            p += n
        while j < len(sites) and sites[j] < p:
            j += 1
    return out


# ---------------------------------------------------------------- het
def het(a):
    loci_p, reg_p, f1, f2 = a
    subs = defaultdict(Counter)
    refb = {}
    als = load_paf(f1, "dep1") + load_paf(f2, "dep2")
    for q, ql, qs, qe, st, t, ts, te, cs, asm in als:
        p = ts
        for m in tok.findall(cs):
            if m[0] == ":":
                p += int(m[1:])
            elif m[0] == "*":
                subs[(asm, t, p)][m[2].upper()] += 1
                refb[(asm, t, p)] = m[1].upper()
                p += 1
            elif m[0] == "-":
                p += len(m) - 1
    cand = defaultdict(list)
    for (asm, t, p), d in subs.items():
        if max(d.values()) >= MIN_ALT:
            cand[(asm, t)].append(p)
    for k in cand:
        cand[k].sort()
    depth = Counter()
    for q, ql, qs, qe, st, t, ts, te, cs, asm in als:
        v = cand.get((asm, t))
        if not v:
            continue
        for p in v[bisect.bisect_left(v, ts):bisect.bisect_left(v, te)]:
            depth[(asm, t, p)] += 1
    print("asm\tchrom\tpos\tref\talt\tdepth\taltn")
    for (asm, t), v in sorted(cand.items()):
        for p in v:
            d = subs[(asm, t, p)]
            alt, n = d.most_common(1)[0]
            dp = depth[(asm, t, p)]
            if dp >= MIN_DEPTH and LO <= n / dp <= HI:
                print(f"{asm}\t{t}\t{p}\t{refb[(asm, t, p)]}\t{alt}\t{dp}\t{n}")


def load_loci(p):
    rows = [l.rstrip("\n").split("\t") for l in open(p)][1:]
    return [(r[0], r[1], r[2], r[3], int(r[4]), int(r[5]), r[6]) for r in rows]


def load_het(p):
    h = defaultdict(dict)
    for l in list(open(p))[1:]:
        asm, c, pos, ref, alt, dp, n = l.rstrip("\n").split("\t")
        h[(asm, c)][int(pos)] = (ref, alt)
    return h


# ---------------------------------------------------------------- phase (test A)
def phase(a):
    loci_p, het_p, f1, f2, outp = a
    loci = load_loci(loci_p)
    het_ = load_het(het_p)
    hpos = {k: sorted(v) for k, v in het_.items()}
    als = load_paf(f1, "dep1") + load_paf(f2, "dep2")
    by_t = defaultdict(list)
    for x in als:
        by_t[(x[9], x[5])].append(x)
    for k in by_t:
        by_t[k].sort(key=lambda x: x[6])
    jrows, lrows = [], []
    for lid, cat, asm, c, ps, pe, cause in loci:
        hp = hpos.get((asm, c), [])
        H = het_.get((asm, c), {})
        piece_het = [x for x in hp if ps <= x < pe and (x - ps < F or pe - x <= F)]
        piece_scan = min(pe - ps, 2 * F)
        flank_het = [x for x in hp if (ps - F <= x < ps) or (pe <= x < pe + F)]
        verdicts = []
        for side, J in (("L", ps), ("R", pe)):
            fl = [x for x in hp if (J - F <= x < J if side == "L" else J <= x < J + F)]
            near = sum(1 for x in fl if abs(x - J) <= NEAR)
            reads = defaultdict(lambda: {"ref": 0, "alt": 0, "cross": False, "stop": None, "dist": None})
            alist = by_t.get((asm, c), [])
            for q, ql, qs, qe, st, t, ts, te, cs, _ in alist:
                if ts > J + F or te < J - F:
                    continue
                cross = ts <= J - J_IN and te >= J + J_IN
                if side == "L":           # approach: flank alignment ending within APPROACH of J (dist = J - end)
                    dist = J - te
                    appr = (not cross) and -J_IN < dist <= APPROACH and ts <= J - J_IN
                    rest = (ql - qe) if st == "+" else qs
                else:
                    dist = ts - J
                    appr = (not cross) and -J_IN < dist <= APPROACH and te >= J + J_IN
                    rest = qs if st == "+" else (ql - qe)
                stop = appr and dist <= STOP_WIN
                if not (cross or appr):
                    continue
                sites = [x for x in fl[bisect.bisect_left(fl, ts):bisect.bisect_left(fl, te)]]
                al = walk(ts, cs, sites)
                r = reads[q]
                for x, b in al.items():
                    if b == "=":
                        r["ref"] += 1
                    elif b == H[x][1]:
                        r["alt"] += 1
                if cross:
                    r["cross"] = True
                else:
                    if r["dist"] is None or dist < r["dist"]:
                        r["dist"] = dist
                    if stop and r["stop"] is None:
                        r["stop"] = "clip" if rest >= 1000 else "readend"
            cnt = Counter()
            closest = {}
            for q, r in reads.items():
                n = r["ref"] + r["alt"]
                g = "unph" if n == 0 else ("ref" if r["ref"] / n >= 0.8 else ("alt" if r["alt"] / n >= 0.8 else "mix"))
                if r["cross"]:
                    cnt[(g, "cross")] += 1
                    continue
                cnt[(g, "approach")] += 1
                closest[g] = min(closest.get(g, 10 ** 9), r["dist"])
                if r["stop"]:
                    cnt[(g, "stop_" + r["stop"])] += 1
            rc_, ac_ = cnt[("ref", "cross")], cnt[("alt", "cross")]
            if min(rc_, ac_) >= 2:
                v = "both"
            elif min(rc_, ac_) == 0 and max(rc_, ac_) >= 4:
                v = "one"
            elif max(rc_, ac_) + min(rc_, ac_) >= 2:
                v = "weak"
            else:
                v = "untestable"
            verdicts.append(v)
            jrows.append([lid, cat, asm, c, side, J, len(fl), near, rc_, ac_, cnt[("unph", "cross")], cnt[("mix", "cross")],
                          cnt[("ref", "stop_clip")], cnt[("alt", "stop_clip")], cnt[("ref", "stop_readend")],
                          cnt[("alt", "stop_readend")], cnt[("unph", "stop_clip")], cnt[("ref", "approach")],
                          cnt[("alt", "approach")], closest.get("ref", "-"), closest.get("alt", "-"), v])
        if "both" in verdicts and "one" in verdicts:
            lv = "conflict"
        elif "both" in verdicts:
            lv = "two_copy"
        elif "one" in verdicts:
            lv = "one_copy"
        elif "weak" in verdicts:
            lv = "weak"
        else:
            lv = "untestable"
        pk = piece_scan / 1000
        fk = 2 * F / 1000
        lrows.append([lid, cat, asm, c, ps, pe, pe - ps, cause, len(flank_het), f"{len(flank_het) / fk:.3f}",
                      len(piece_het), f"{len(piece_het) / pk:.3f}" if pk else "-", verdicts[0], verdicts[1], lv])
    with open(outp + ".junctions.tsv", "w") as fh:
        fh.write("locus\tcat\tasm\tchrom\tside\tjunction\tflank_het\tflank_het_10kb\tref_cross\talt_cross\tunph_cross\t"
                 "mix_cross\tref_stop_clip\talt_stop_clip\tref_stop_readend\talt_stop_readend\tunph_stop_clip\t"
                 "ref_approach20k\talt_approach20k\tref_closest_end\talt_closest_end\tverdict\n")
        for r in jrows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(outp + ".loci.tsv", "w") as fh:
        fh.write("locus\tcat\tasm\tchrom\tps\tpe\tlen\tcause\tflank_het\tflank_het_per_kb\tpiece_het\tpiece_het_per_kb\t"
                 "vL\tvR\tverdict\n")
        for r in lrows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(outp + ".summary.txt", "w") as fh:
        for cat in ("cat1", "cat2", "ctl_rand", "ctl_het"):
            rs = [r for r in lrows if r[1] == cat]
            if not rs:
                continue
            fh.write(f"{cat}\tn={len(rs)}\t" + " ".join(f"{k}:{v}" for k, v in sorted(Counter(r[-1] for r in rs).items()))
                     + f"\tpiece_het>0:{sum(1 for r in rs if r[10] > 0)}"
                     + f"\tmedian_flank_het_per_kb:{sorted(float(r[9]) for r in rs)[len(rs) // 2]:.3f}\n")


# ---------------------------------------------------------------- Hi-C (test B)
def hicanc(a):
    loci_p, het_p, regfa, out = a
    loci = load_loci(loci_p)
    het_ = load_het(het_p)
    regs = defaultdict(list)
    for n, s in read_fa(regfa).items():
        asm, rest = n.split("|")
        c, iv = rest.rsplit(":", 1)
        b, e = map(int, iv.split("-"))
        regs[(asm, c)].append((b, e, s))
    want = set()
    for lid, cat, asm, c, ps, pe, cause in loci:
        for x in het_.get((asm, c), {}):
            if (ps - F <= x < ps) or (pe <= x < pe + F):
                want.add((asm, c, x))
    with open(f"{out}/hic_sites.tsv", "w") as fs, open(f"{out}/hic_cand.fa", "w") as fc:
        n = 0
        for asm, c, x in sorted(want):
            hp = het_[(asm, c)]
            if any(abs(y - x) < K_HIC and y != x for y in range(x - K_HIC, x + K_HIC) if y in hp):
                continue
            seq = next(((s, b) for b, e, s in regs[(asm, c)] if b + K_HIC <= x < e - K_HIC), None)
            if seq is None:
                continue
            s, b = seq
            i = x - b
            ref, alt = hp[x]
            if s[i] != ref:
                continue
            fs.write(f"{asm}\t{c}\t{x}\t{ref}\t{alt}\n")
            for allele in (ref, alt):
                w = s[i - K_HIC + 1:i] + allele + s[i + 1:i + K_HIC]
                for o in range(K_HIC):
                    k = w[o:o + K_HIC]
                    if "N" not in k:
                        n += 1
                        fc.write(f">{asm}|{c}|{x}|{allele}|{o}\n{k}\n")


def hic(a):
    loci_p, sites_p, paf1, paf2, outp = a
    loci = load_loci(loci_p)
    sites = defaultdict(dict)
    for l in open(sites_p):
        asm, c, x, ref, alt = l.split()
        sites[(asm, c)][int(x)] = (ref, alt)
    spos = {k: sorted(v) for k, v in sites.items()}
    # per asm: read -> mate -> (chrom, ts, te, mapq, alleles)
    res = defaultdict(lambda: defaultdict(list))
    for asm, paf in (("dep1", paf1), ("dep2", paf2)):
        with opn(paf) as fh:
            for l in fh:
                f = l.rstrip("\n").split("\t")
                if "tp:A:P" not in f[12:]:
                    continue
                q, mate = f[0].rsplit("_", 1)
                t, ts, te, mq = f[5], int(f[7]), int(f[8]), int(f[11])
                cs = next(x for x in f[12:] if x.startswith("cs:Z:"))[5:]
                sp = spos.get((asm, t), [])
                al = walk(ts, cs, sp[bisect.bisect_left(sp, ts):bisect.bisect_left(sp, te)]) if sp else {}
                res[asm][q].append((mate, t, ts, te, mq, al))
    # locus lookup per (asm, chrom)
    lk = defaultdict(list)
    for L in loci:
        lk[(L[2], L[3])].append(L)
    cnt = defaultdict(Counter)   # (locus, side, allele) -> class counts
    for asm in res:
        for q, ms in res[asm].items():
            mates = {}
            for m in ms:
                if m[0] not in mates or m[4] > mates[m[0]][4]:
                    mates[m[0]] = m
            if len(mates) != 2:
                continue
            for a_m, b_m in ((mates["1"], mates["2"]), (mates["2"], mates["1"])):
                mate, t, ts, te, mq, al = a_m
                _, t2, ts2, te2, mq2, _ = b_m
                if mq < 20 or mq2 < 20 or t2 != t or not al:
                    continue
                qpos = (ts2 + te2) // 2
                for x, b in al.items():
                    ref, alt = sites[(asm, t)][x]
                    allele = "ref" if b == "=" else ("alt" if b == alt else None)
                    if allele is None:
                        continue
                    for lid, cat, _, c, ps, pe, cause in lk[(asm, t)]:
                        if ps - F <= x < ps:
                            side = "L"
                        elif pe <= x < pe + F:
                            side = "R"
                        else:
                            continue
                        if abs(qpos - x) < 1000:
                            k = "self"
                        elif ps <= qpos < pe:
                            k = "piece"
                        elif (side == "L" and qpos < ps and qpos >= ps - F) or (side == "R" and qpos >= pe and qpos < pe + F):
                            k = "same"
                        elif (side == "L" and pe <= qpos < pe + F) or (side == "R" and ps - F <= qpos < ps):
                            k = "other"
                        else:
                            k = "far"
                        cnt[(lid, side, allele)][k] += 1
    rows = []
    for lid, cat, asm, c, ps, pe, cause in loci:
        sv = []
        for side in ("L", "R"):
            r, al_ = cnt[(lid, side, "ref")], cnt[(lid, side, "alt")]
            pr, sr, pa, sa = r["piece"], r["same"], al_["piece"], al_["same"]
            fr = pr / (pr + sr) if pr + sr else float("nan")
            fa = pa / (pa + sa) if pa + sa else float("nan")
            p = fisher(pr, sr, pa, sa)
            if pr + sr < 20 or pa + sa < 20 or pr + pa < 10:
                v = "untestable"
            else:
                lo, hi = sorted((fr, fa))
                ratio = lo / hi if hi > 0 else float("nan")
                if p < 0.001 and ratio < 0.3:
                    v = "one"
                elif ratio >= 0.5:
                    v = "both"
                else:
                    v = "unclear"
            sv.append(v)
            rows.append([lid, cat, asm, side, pe - ps, pr, sr, r["other"], pa, sa, al_["other"],
                         f"{fr:.3f}", f"{fa:.3f}", f"{p:.2g}", v])
        if "both" in sv and "one" in sv:
            lv = "conflict"
        elif "both" in sv:
            lv = "two_copy"
        elif "one" in sv:
            lv = "one_copy"
        elif "unclear" in sv:
            lv = "unclear"
        else:
            lv = "untestable"
        for r in rows[-2:]:
            r.append(lv)
    with open(outp + ".tsv", "w") as fh:
        fh.write("locus\tcat\tasm\tside\tlen\tref_piece\tref_same\tref_other\talt_piece\talt_same\talt_other\t"
                 "ref_frac\talt_frac\tfisher_p\tside_verdict\tlocus_verdict\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(outp + ".summary.txt", "w") as fh:
        for cat in ("cat1", "cat2", "ctl_rand", "ctl_het"):
            ls = {r[0]: r[-1] for r in rows if r[1] == cat}
            if ls:
                fh.write(f"{cat}\tn={len(ls)}\t" + " ".join(f"{k}:{v}" for k, v in sorted(Counter(ls.values()).items())) + "\n")
            # pooled piece fractions per category (allele with the lower fraction per side vs the higher)
            lo, hi_ = [0, 0], [0, 0]
            for r in rows:
                if r[1] != cat or r[-2] == "untestable":
                    continue
                a_ = (r[5], r[6]); b_ = (r[8], r[9])
                if a_[0] / max(1, sum(a_)) > b_[0] / max(1, sum(b_)):
                    a_, b_ = b_, a_
                lo = [lo[0] + a_[0], lo[1] + a_[1]]; hi_ = [hi_[0] + b_[0], hi_[1] + b_[1]]
            if sum(lo) and sum(hi_):
                fh.write(f"{cat}\tpooled testable sides: lower-allele piece frac {lo[0] / sum(lo):.3f} ({lo[0]}/{sum(lo)}),"
                         f" higher-allele {hi_[0] / sum(hi_):.3f} ({hi_[0]}/{sum(hi_)})\n")


def fisher(a, b, c, d):
    """two-sided Fisher exact p for [[a,b],[c,d]]"""
    n = a + b + c + d
    r1, c1 = a + b, a + c
    if n == 0:
        return 1.0
    lg = lambda k: math.lgamma(k + 1)
    base = lg(r1) + lg(n - r1) + lg(c1) + lg(n - c1) - lg(n)
    pr = lambda x: math.exp(base - lg(x) - lg(r1 - x) - lg(c1 - x) - lg(n - r1 - c1 + x))
    p0 = pr(a)
    lo, hi = max(0, r1 + c1 - n), min(r1, c1)
    return min(1.0, sum(pr(x) for x in range(lo, hi + 1) if pr(x) <= p0 * (1 + 1e-7)))


# ---------------------------------------------------------------- betulina (test C)
def bet(a):
    """PAF of query '<locus>|<qps>|<qpe>' (piece +-10 kb, piece at qps..qpe) vs betulina primary (asm20, -c).
    Per locus: best target = most aligned query bases; 'carries' when that target's same-strand blocks cover >= 80 %
    of the piece and >= 50 % of each 10 kb flank, collinear (target order follows query order)."""
    loci_p, paf = a
    al = defaultdict(list)
    for l in opn(paf):
        f = l.split("\t")
        al[f[0]].append((int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])))
    print("locus\tcat\tlen\ttarget\tstrand\tpiece_cov\tflankL_cov\tflankR_cov\tcollinear\tident\tverdict")
    for lid, cat, asm, c, ps, pe, cause in load_loci(loci_p):
        qn = next((q for q in al if q.split("|")[0] == lid), None)
        if qn is None:
            print(f"{lid}\t{cat}\t{pe - ps}\t-\t-\t0\t0\t0\t-\t-\tno_hit")
            continue
        _, qps, qpe = qn.split("|"); qps, qpe = int(qps), int(qpe)
        best = Counter()
        for x in al[qn]:
            best[(x[3], x[2])] += x[1] - x[0]
        (tg, st), _ = best.most_common(1)[0]
        bl = sorted(x for x in al[qn] if x[3] == tg and x[2] == st)
        cov = lambda s, e: sum(max(0, min(e, x[1]) - max(s, x[0])) for x in merge_q(bl)) / max(1, e - s)
        pc, lc, rc_ = cov(qps, qpe), cov(0, qps), cov(qpe, qpe + qps)
        tpos = [x[4] for x in bl]
        coll = tpos == sorted(tpos) if st == "+" else tpos == sorted(tpos, reverse=True)
        ident = sum(x[6] for x in bl) / max(1, sum(x[7] for x in bl))
        v = "carries" if pc >= 0.8 and lc >= 0.5 and rc_ >= 0.5 and coll else ("partial" if pc >= 0.5 else "absent")
        print(f"{lid}\t{cat}\t{pe - ps}\t{tg}\t{st}\t{pc:.2f}\t{lc:.2f}\t{rc_:.2f}\t{coll}\t{ident:.3f}\t{v}")


def merge_q(bl):
    return [tuple(x) for x in merge([(x[0], x[1]) for x in bl])]


def betq(a):
    """write query FASTA for bet: piece +-10 kb from the gene-search assemblies"""
    loci_p, gs1, gs2, maxlen = a[0], a[1], a[2], int(a[3])
    gs = {"dep1": read_fa(gs1), "dep2": read_fa(gs2)}
    for lid, cat, asm, c, ps, pe, cause in load_loci(loci_p):
        if pe - ps > maxlen:
            continue
        s = max(0, ps - 10000)
        seq = gs[asm][c][s:pe + 10000]
        print(f">{lid}|{ps - s}|{pe - s}\n{seq}")


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    {"prep": prep, "anchors": anchors, "filter": filt, "het": het, "phase": phase, "hicanc": hicanc, "hic": hic,
     "bet": bet, "betq": betq}[cmd](args)
