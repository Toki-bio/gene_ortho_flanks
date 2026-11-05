#!/usr/bin/env python3
"""Phase test of fills. Is a fill this haplotype's own sequence, or the other haplotype's allele?
For each test: target window T (left flank + N run + right flank), other-haplotype window O, one or more fixed windows
(left flank + fill + right flank). Each read is mapped to each window ALONE.
Assignment: difference events on T and on O, counted over the same read bases, flanks only (never inside the hole).
Crossing: read aligned on T within 500 bp on both sides of the hole, same strand.
Fill check: differences inside the fill on the fixed window.

usage: phase_fill.py controls   -> TMED4 + 20 ordinary sites, 2 kb hidden behind 100 N, fills = own and other allele
"""
import os, random, re, statistics, subprocess, sys, collections
G = "/home/toki/gof"; E = "/home/toki/miniforge3/envs"
S = f"{E}/sinederella/bin/samtools"; M = f"{E}/hifiasm/bin/minimap2"
T1 = "/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna"
T2 = "/home/toki/Sicista2026/demog/genomes/tri_hap2.fna"
HC = f"{G}/asm2/hicanu/tri.contigs.fasta"; HB = f"{G}/asm2/hicanu/tri.contigs.bam"
FISHED = f"{G}/repair/fished.fq.gz"
W = f"{G}/phase_test"; TH = 16
FL = 25000          # flank
BIG = 50            # indels from this size are not counted as assignment events
SHARED = 2000       # minimum shared read bases for an assignment
MARGIN = 3          # minimum difference in events
OWN = 0.005         # own-haplotype rate above this = read from elsewhere
EDGE = 500; EDGE_MIN = 400
CS = re.compile(r":(\d+)|\*([a-z])([a-z])|\+([a-z]+)|-([a-z]+)")


def sh(c):
    return subprocess.run(c, shell=True, capture_output=True, text=True).stdout


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))


def get(fa, c, a, b):
    """0-based half-open"""
    return "".join(sh(f"{S} faidx {fa} {c}:{a+1}-{b}").split("\n")[1:]).upper()


def wfa(p, d):
    with open(p, "w") as o:
        for k, v in d.items():
            o.write(f">{k}\n{v}\n")


class Aln:
    """All alignments of each read on one window."""

    def __init__(self, paf, hole=None):
        self.r = collections.defaultdict(list)
        for l in sh(f"cat {paf}").splitlines():
            f = l.split("\t")
            cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
            if cs is None:
                continue
            self.r[f[0]].append(dict(L=int(f[1]), qs=int(f[2]), qe=int(f[3]), st=f[4], rs=int(f[7]), re=int(f[8]), cs=cs))
        self.hole = hole

    def walk(self, a):
        """yield (kind, ref, read_low, n): kind m=matched run, x=mismatch, i=insertion in read, d=deletion in read.
        read_low = lowest read coordinate of the run."""
        r = a["rs"]; o = 0; fw = a["st"] == "+"; span = a["qe"] - a["qs"]
        for m in CS.finditer(a["cs"]):
            if m.group(1):
                n = int(m.group(1))
                yield "m", r, (a["qs"] + o if fw else a["qe"] - o - n), n
                r += n; o += n
            elif m.group(2):
                yield "x", r, (a["qs"] + o if fw else a["qe"] - o - 1), 1
                r += 1; o += 1
            elif m.group(4):
                n = len(m.group(4))
                yield "i", r, (a["qs"] + o if fw else a["qe"] - o - n), n
                o += n
            else:
                n = len(m.group(5))
                yield "d", r, (a["qs"] + o if fw else a["qe"] - o), n
                r += n

    def cover(self, read):
        """read bases aligned outside the hole, plus events there. Returns bytearray, events, per-side edge cover and strand."""
        al = self.r.get(read, [])
        if not al:
            return None
        cov = bytearray(al[0]["L"] + 1); ev = []; left = {"+": 0, "-": 0}; right = {"+": 0, "-": 0}
        g0, g1 = self.hole if self.hole else (None, None)
        for a in al:
            for k, r, q, n in self.walk(a):
                if k in "mx":
                    for j in range(n):
                        p = r + j
                        rq = q + j if a["st"] == "+" else q + n - 1 - j
                        if g0 is not None and g0 <= p < g1:
                            continue
                        cov[rq] = 1
                        if g0 is not None:
                            if g0 - EDGE <= p < g0:
                                left[a["st"]] += 1
                            elif g1 <= p < g1 + EDGE:
                                right[a["st"]] += 1
                    if k == "x" and not (g0 is not None and g0 <= r < g1):
                        ev.append(q)
                elif n < BIG and not (g0 is not None and g0 - 1 <= r <= g1):
                    ev.append(max(q - 1, 0))
        return cov, ev, left, right

    def in_fill(self, read, f0, f1):
        """differences inside [f0,f1) on a fixed window."""
        out = dict(cov=0, x=0, ind=0, indbp=0, big=[], spans=False)
        for a in self.r.get(read, []):
            if a["rs"] <= f0 - EDGE_MIN and a["re"] >= f1 + EDGE_MIN:
                out["spans"] = True
            for k, r, q, n in self.walk(a):
                if k == "m":
                    out["cov"] += max(0, min(r + n, f1) - max(r, f0))
                elif f0 <= r < f1:
                    if k == "x":
                        out["x"] += 1; out["cov"] += 1
                    else:
                        out["ind"] += 1; out["indbp"] += n
                        if n >= BIG:
                            out["big"].append(("+" if k == "i" else "-") + str(n))
        return out


def assign(t, o):
    """t, o = cover() results. Returns call, shared, dT, dO"""
    if t is None or o is None:
        return "one_side", 0, 0, 0
    ct, et = t[0], t[1]; co, eo = o[0], o[1]
    shared = sum(1 for x, y in zip(ct, co) if x and y)
    dt = sum(1 for p in et if ct[p] and co[p]); do = sum(1 for p in eo if ct[p] and co[p])
    if shared < SHARED:
        return "short", shared, dt, do
    if min(dt, do) / shared > OWN:
        return "foreign", shared, dt, do
    if do - dt >= MARGIN and dt * 2 <= do:
        return "T", shared, dt, do
    if dt - do >= MARGIN and do * 2 <= dt:
        return "O", shared, dt, do
    return "tie", shared, dt, do


def paf(ref, fq, out):
    sh(f"{M} -x map-hifi -c --cs=short --secondary=no -t {TH} {ref} {fq} 2>/dev/null > {out}")


def run_test(t, fq):
    """t: id, T, hole, O, fills {name: seq}"""
    d = f"{W}/{t['id']}"; os.makedirs(d, exist_ok=True)
    wfa(f"{d}/T.fa", {"T": t["T"]}); wfa(f"{d}/O.fa", {"O": t["O"]})
    paf(f"{d}/T.fa", fq, f"{d}/T.paf"); paf(f"{d}/O.fa", fq, f"{d}/O.paf")
    g0, g1 = t["hole"]
    aT = Aln(f"{d}/T.paf", (g0, g1)); aO = Aln(f"{d}/O.paf")
    fx = {}
    for n, s in t["fills"].items():
        wfa(f"{d}/F_{n}.fa", {"F": t["T"][:g0] + s + t["T"][g1:]})
        paf(f"{d}/F_{n}.fa", fq, f"{d}/F_{n}.paf")
        fx[n] = (Aln(f"{d}/F_{n}.paf"), g0, g0 + len(s))
    calls = collections.Counter(); cross = collections.Counter(); rows = []
    res = {n: collections.defaultdict(list) for n in fx}
    # flank difference between the two windows, as seen by the reads: reported through dT and dO
    for rd in aT.r:
        ct = aT.cover(rd); co = aO.cover(rd)
        call, shared, dt, do = assign(ct, co)
        calls[call] += 1
        st = "+" if ct[2]["+"] + ct[3]["+"] >= ct[2]["-"] + ct[3]["-"] else "-"
        x = ct[2][st] >= EDGE_MIN and ct[3][st] >= EDGE_MIN
        if x:
            cross[call] += 1
        extra = []
        if x:
            for n, (a, f0, f1) in fx.items():
                r = a.in_fill(rd, f0, f1)
                res[n][call].append(r)
                extra.append(f"{n}:cov={r['cov']},x={r['x']},ind={r['ind']},indbp={r['indbp']},big={'/'.join(r['big'])},spans={int(r['spans'])}")
        rows.append([t["id"], rd, call, shared, dt, do, int(x), ";".join(extra)])
    return calls, cross, res, rows


def per_kb(rs):
    v = [1000 * (r["x"] + r["ind"]) / max(r["cov"], 1) for r in rs if r["spans"]]
    return (f"{statistics.median(v):.2f}", len(v)) if v else ("NA", 0)


def controls():
    tests = []
    # TMED4: real gap, real fill
    g0, g1 = 96007474, 96007574; c1 = "OZ418347.1"; c2 = "OZ418087.1"; s0, s1 = 95802104, 95804435
    Tw = get(T1, c1, g0 - FL, g1 + FL)
    Ow = get(T2, c2, s0 - FL, s1 + FL)
    assert set(Tw[FL:FL + 100]) == {"N"}, "TMED4 hole is not N"
    tests.append(dict(id="TMED4", T=Tw, hole=(FL, FL + 100), O=Ow, fills={"real": get(T2, c2, s0, s1)},
                      tloc=(c1, g0 - FL, g1 + FL), oloc=(c2, s0 - FL, s1 + FL), strand="+"))
    # ordinary sites
    random.seed(1105)
    fai = [l.split("\t") for l in open(T1 + ".fai")]
    chrom = [(f[0], int(f[1])) for f in fai if int(f[1]) > 20_000_000]
    HOLE = 2000; cand = {}
    while len(cand) < 40:
        c, n = random.choices(chrom, weights=[x[1] for x in chrom])[0]
        p = random.randrange(1_000_000, n - 1_000_000)
        s = get(T1, c, p - FL, p + HOLE + FL)
        if "N" not in s:
            cand[f"ctl_{c}_{p}"] = (c, p, s)
    wfa(f"{W}/cand.fa", {k: v[2] for k, v in cand.items()})
    o = sh(f"{M} -x asm20 -c --cs=short -r 2000,120000 -z 1000000 --secondary=no -t {TH} {T2} {W}/cand.fa 2>/dev/null")
    a = {}
    for l in o.splitlines():
        f = l.split("\t")
        if int(f[3]) - int(f[2]) >= 40000 and f[0] not in a and int(f[2]) < FL - 1000 and int(f[3]) > FL + HOLE + 1000:
            a[f[0]] = f
    # hidden 2 kb, and its position on hap2 by walking the alignment
    for k, f in a.items():
        if len(tests) > 20:
            break
        c, p, s = cand[k]
        st = f[4]; cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
        qs, qe, rs = int(f[2]), int(f[3]), int(f[7])
        q = qs if st == "+" else qe; r = rs; want = {FL, FL + HOLE}; got = {}
        for m in CS.finditer(cs):
            if m.group(1):
                n = int(m.group(1))
            elif m.group(2):
                n = 1
            elif m.group(4):
                q += len(m.group(4)) if st == "+" else -len(m.group(4)); continue
            else:
                r += len(m.group(5)); continue
            for w in list(want):
                off = (w - q) if st == "+" else (q - w)
                if 0 <= off < n:
                    got[w] = r + off; want.discard(w)
            q += n if st == "+" else -n; r += n
        if len(got) < 2:
            continue
        lo, hi = sorted((got[FL], got[FL + HOLE]))
        Ow = get(T2, f[5], lo - FL, hi + FL)
        if len(Ow) < 2 * FL or "N" in Ow:
            continue
        other = get(T2, f[5], lo, hi)
        if st == "-":
            other = rc(other); Ow = rc(Ow)
        own = s[FL:FL + HOLE]
        tests.append(dict(id=k, T=s[:FL] + "N" * 100 + s[FL + HOLE:], hole=(FL, FL + 100), O=Ow, fills={"own": own, "other": other},
                          tloc=(c, p - FL, p + HOLE + FL), oloc=(f[5], lo - FL, hi + FL), strand=st))
    print("tests", len(tests), flush=True)
    # reads: the contigs under both windows, once
    q = {}
    for t in tests:
        q[t["id"] + "|T"] = t["T"].replace("N", ""); q[t["id"] + "|O"] = t["O"]
    wfa(f"{W}/wins.fa", q)
    o = sh(f"{M} -x asm10 -c --secondary=no -t {TH} {HC} {W}/wins.fa 2>/dev/null")
    reg = collections.defaultdict(set)
    for l in o.splitlines():
        f = l.split("\t")
        if int(f[9]) >= 2000:
            reg[f[0].split("|")[0]].add(f"{f[5]}:{int(f[7])+1}-{f[8]}")
    out = open(f"{W}/controls.tsv", "w"); rr = open(f"{W}/controls.reads.tsv", "w")
    print("id\tT\tO\tstrand\treads\tT_reads\tO_reads\ttie\tforeign\tshort\tone_side\tT_cross\tO_cross\ttie_cross\tfill\tT_cross_diff_per_kb(n)\tO_cross_diff_per_kb(n)", file=out)
    for t in tests:
        d = f"{W}/{t['id']}"; os.makedirs(d, exist_ok=True)
        fq = f"{d}/reads.fq"
        sh(f"{S} view -b -F 0x900 {HB} {' '.join(sorted(reg[t['id']]))} | {S} fastq - 2>/dev/null > {d}/hc.fq")
        src = f"cat {d}/hc.fq; zcat {FISHED}" if t["id"] == "TMED4" else f"cat {d}/hc.fq"
        seen = set()
        p = subprocess.Popen(src, shell=True, stdout=subprocess.PIPE, text=True)
        with open(fq, "w") as w:
            while True:
                h = p.stdout.readline()
                if not h:
                    break
                b = [p.stdout.readline() for _ in range(3)]
                n = h.split()[0]
                if n not in seen:
                    seen.add(n); w.write(h); w.writelines(b)
        calls, cross, res, rows = run_test(t, fq)
        for r in rows:
            print(*r, sep="\t", file=rr)
        tot = sum(calls.values())
        for n in t["fills"]:
            a = per_kb(res[n]["T"]); b = per_kb(res[n]["O"])
            print(t["id"], "%s:%d-%d" % t["tloc"], "%s:%d-%d" % t["oloc"], t["strand"], tot, calls["T"], calls["O"], calls["tie"], calls["foreign"], calls["short"], calls["one_side"],
                  cross["T"], cross["O"], cross["tie"], n, f"{a[0]}({a[1]})", f"{b[0]}({b[1]})", sep="\t", file=out)
        out.flush(); rr.flush()
        for f in os.listdir(d):
            if f.endswith((".fq", ".paf")) and t["id"] != "TMED4":
                os.remove(f"{d}/{f}")
        print(t["id"], "done", dict(calls), "cross", dict(cross), flush=True)
    print("finished", flush=True)


if __name__ == "__main__":
    os.makedirs(W, exist_ok=True)
    if sys.argv[1] == "controls":
        controls()
