#!/usr/bin/env python3
"""Phase test of fills. Is a fill this haplotype's own sequence, or the other haplotype's allele?
For each test: target window T (left flank + N run + right flank), other-haplotype window O, one or more fixed windows
(left flank + fill + right flank). Each read is mapped to each window ALONE.
Assignment: difference events on T and on O, counted over the same read bases, flanks only (never inside the hole).
Crossing: read aligned on T within 500 bp on both sides of the hole, same strand.
Fill check: differences inside the fill on the fixed window.

usage: phase_fill.py controls   -> TMED4 + 20 ordinary sites, 2 kb hidden behind 100 N, fills = own and other allele
       phase_fill.py pool       -> all HiFi reads mapped once to each deposited haplotype alone, reads in the fill windows kept
       phase_fill.py all        -> every accepted fill; fixed window = slice of repair/depN.genesearch.fa; stops if TMED4 is not wrong_allele
Verdict uses substitutions inside the fill, target reads that cross and span it: median 0 own_or_identical, >=2 wrong_allele, else undecided.
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
                yield ("n" if "n" in m.group(2) + m.group(3) else "x"), r, (a["qs"] + o if fw else a["qe"] - o - 1), 1
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
                if k == "n":
                    continue
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
                if k == "n":
                    continue
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
    for n, s in t.get("fills", {}).items():
        wfa(f"{d}/F_{n}.fa", {"F": t["T"][:g0] + s + t["T"][g1:]})
        paf(f"{d}/F_{n}.fa", fq, f"{d}/F_{n}.paf")
        fx[n] = (Aln(f"{d}/F_{n}.paf"), g0, g0 + len(s))
    for n, (s, f0, f1) in t.get("fixed", {}).items():
        wfa(f"{d}/F_{n}.fa", {"F": s})
        paf(f"{d}/F_{n}.fa", fq, f"{d}/F_{n}.paf")
        fx[n] = (Aln(f"{d}/F_{n}.paf"), f0, f1)
    calls = collections.Counter(); cross = collections.Counter(); rows = []
    res = {n: collections.defaultdict(list) for n in fx}
    enter = {n: collections.defaultdict(list) for n in fx}
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
        for n, (a, f0, f1) in fx.items():
            r = a.in_fill(rd, f0, f1)
            if not x:
                if r["cov"] >= EDGE:
                    enter[n][call].append(r)
                    extra.append(f"{n}:enter,cov={r['cov']},x={r['x']},ind={r['ind']},indbp={r['indbp']}")
            else:
                res[n][call].append(r)
                extra.append(f"{n}:cov={r['cov']},x={r['x']},ind={r['ind']},indbp={r['indbp']},big={'/'.join(r['big'])},spans={int(r['spans'])}")
        rows.append([t["id"], rd, call, shared, dt, do, int(x), ";".join(extra)])
    return calls, cross, res, rows, enter


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
        calls, cross, res, rows, _ = run_test(t, fq)
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


# ---------- all accepted fills ----------
GEN = {"dep1": T1, "dep2": T2}
READS = f"{G}/reads/ERR16751478.fastq.gz"
HEAD = ("assembly gap chrom start end cause action ins_len source_name source strand reads T_reads O_reads tie foreign short one_side "
        "T_cross O_cross tie_cross T_span T_x T_x_median T_indel_median T_big O_x_median tie_x_median T_enter T_enter_x_per_kb O_enter_x_per_kb verdict").split()
POOL = {}


def verdict(calls, cross, rs):
    if cross["T"] == 0:
        if cross["O"] == 0 and cross["tie"] > 0:
            return "flanks_tie"
        if calls["T"] + calls["O"] == 0:
            return "flanks_tie" if calls["tie"] else "no_reads"
        return "no_target_cross"
    sp = [r for r in rs if r["spans"]]
    if not sp:
        return "wrong_allele_structure"
    m = statistics.median(r["x"] for r in sp)
    return "own_or_identical" if m == 0 else "wrong_allele" if m >= 2 else "undecided"


def med(rs, k="x"):
    v = [r[k] for r in rs if r["spans"]]
    return statistics.median(v) if v else "NA"


def ekb(rs):
    c = sum(r["cov"] for r in rs)
    return f"{1000 * sum(r['x'] for r in rs) / c:.2f}" if c else "NA"


def one(t):
    d = f"{W}/run/{t['id']}"; os.makedirs(d, exist_ok=True)
    try:
        t = dict(t, id=f"run/{t['id']}")
        fq = f"{d}/reads.fq"
        with open(fq, "w") as w:
            for n in t["names"]:
                if n in POOL:
                    w.write(POOL[n])
        calls, cross, res, rows, enter = run_test(t, fq)
        rs = res["built"]["T"]; sp = [r for r in rs if r["spans"]]
        big = sorted({b for r in sp for b in r["big"]})
        row = t["meta"] + [sum(calls.values()), calls["T"], calls["O"], calls["tie"], calls["foreign"], calls["short"], calls["one_side"],
                           cross["T"], cross["O"], cross["tie"], len(sp), ",".join(str(r["x"]) for r in sp), med(rs),
                           statistics.median(r["ind"] for r in sp) if sp else "NA", "/".join(big), med(res["built"]["O"]), med(res["built"]["tie"]),
                           len(enter["built"]["T"]), ekb(enter["built"]["T"]), ekb(enter["built"]["O"]), verdict(calls, cross, rs)]
        return row, [[t["meta"][0], t["meta"][1]] + r[1:] for r in rows if r[2] in ("T", "O") or r[6]]
    except Exception as e:
        return t["meta"] + [f"ERROR {e!r}"], []
    finally:
        sh(f"rm -rf {d}")


def owin(fa, c, o0, o1):
    """other-haplotype window; only its flanks are used for assignment, so a long middle is dropped"""
    a = max(0, o0 - FL)
    if o1 - o0 <= 100000:
        return get(fa, c, a, o1 + FL), [(c, a, o1 + FL)]
    return get(fa, c, a, o0) + "N" * 100 + get(fa, c, o1, o1 + FL), [(c, a, o0), (c, o1, o1 + FL)]


def build_tests():
    r = [l.rstrip("\n").split("\t") for l in open(f"{G}/repair/repairs.tsv")]
    D = [dict(zip(r[0], x)) for x in r[1:]]
    F = [x for x in D if x["accepted"] == "yes" and int(x["ins_len"] or 0) > 0]
    bed = {}
    for a in ("dep1", "dep2"):
        for l in open(f"{G}/repair/{a}.genesearch.bed"):
            f = l.split("\t")
            bed[(a, f[3].split("|")[0])] = (f[0], int(f[1]), int(f[2]))
    tests = []; need = collections.defaultdict(list)
    for x in F:
        a = x["assembly"]; g0, g1 = int(x["start"]), int(x["end"]); c = x["chrom"]
        try:
            l = min(FL, g0)
            Tw = get(GEN[a], c, g0 - l, g1 + FL)
            assert set(Tw[l:l + g1 - g0]) == {"N"}, "hole is not N"
            bc, s, e = bed[(a, x["gap"])]
            assert bc == c and e - s == int(x["ins_len"]), f"bed {bc} {s} {e} {x['ins_len']}"
        except (AssertionError, KeyError) as err:
            print("SKIP", a, x["gap"], repr(err), flush=True)
            continue
        lf = min(FL, s)
        Fw = get(f"{G}/repair/{a}.genesearch.fa", c, s - lf, e + FL)
        src = x["source"].split("|")
        other = "dep2" if a == "dep1" else "dep1"
        t = dict(id=f"{a}_{x['gap']}", T=Tw, hole=(l, l + g1 - g0), fixed={"built": (Fw, lf, lf + e - s)}, O=None, other=other,
                 win={a: [(c, g0 - l, g1 + FL)], other: []},
                 meta=[a, x["gap"], c, g0, g1, x["cause"], x["action"], x["ins_len"], src[1], src[2], src[3]])
        if src[1] == other:
            oc, oab = src[2].rsplit(":", 1); o0, o1 = map(int, oab.split("-"))
            t["O"], t["win"][other] = owin(GEN[other], oc, o0, o1)
        else:
            need[other].append((t, l))
        tests.append(t)
    # fills from hifiasm contigs: find the other deposited haplotype from the flanks
    for other, L in need.items():
        q = {}
        for t, l in L:
            q[t["id"] + "|L"] = t["T"][max(0, l - 15000):l]; q[t["id"] + "|R"] = t["T"][t["hole"][1]:t["hole"][1] + 15000]
        wfa(f"{W}/flank_{other}.fa", {k: v for k, v in q.items() if len(v) >= 1000})
        hit = {}
        for ln in sh(f"{M} -x asm20 --secondary=no -t 32 {GEN[other]} {W}/flank_{other}.fa 2>/dev/null").splitlines():
            f = ln.split("\t")
            if f[0] not in hit or int(f[9]) > int(hit[f[0]][9]):
                hit[f[0]] = f
        for t, l in L:
            a = hit.get(t["id"] + "|L"); b = hit.get(t["id"] + "|R")
            if a and b and a[5] == b[5]:
                x = sorted((int(a[7]), int(a[8]), int(b[7]), int(b[8])))
                t["O"], t["win"][other] = owin(GEN[other], a[5], x[1], x[2])
            elif a or b:
                h = a or b
                t["O"], t["win"][other] = owin(GEN[other], h[5], int(h[7]), int(h[8]))
    return tests


def filt(wfile):
    """stdin PAF of all reads on one genome -> read, test for reads aligned >=1 kb in a window. Keeps no big file."""
    B = 100000; bins = collections.defaultdict(list)
    for l in open(wfile):
        c, a, b, i = l.split()
        a = int(a); b = int(b)
        for k in range(a // B, b // B + 1):
            bins[(c, k)].append((a, b, i))
    seen = set(); n = 0; m = 0
    for l in sys.stdin:
        f = l.split("\t", 11); n += 1
        if int(f[10]) < 1000:
            continue
        c, a, b = f[5], int(f[7]), int(f[8])
        m += 1
        for k in range(a // B, b // B + 1):
            for x, y, i in bins.get((c, k), ()):
                if a < y and b > x and (f[0], i) not in seen:
                    seen.add((f[0], i)); print(f[0], i, sep="\t")
    print("paf lines", n, ">=1 kb blocks", m, "read-test pairs", len(seen), file=sys.stderr)


def pool():
    """Every HiFi read mapped once to each deposited haplotype alone; keep only reads in the windows of the fills."""
    import pickle
    tests = build_tests()
    print("tests", len(tests), "without other window", sum(1 for t in tests if not t["O"]), flush=True)
    pickle.dump(tests, open(f"{W}/tests.pkl", "wb"))
    for g in GEN:
        with open(f"{W}/win_{g}.tsv", "w") as o:
            for t in tests:
                for c, a, b in t["win"].get(g, []):
                    print(c, a, b, t["id"], sep="\t", file=o)
    me = os.path.abspath(__file__)
    ps = [subprocess.Popen(f"set -o pipefail; {M} -x map-hifi --secondary=no -t 40 {GEN[g]} {READS} 2> {W}/map_{g}.log | python3 {me} filter {W}/win_{g}.tsv > {W}/pool_{g}.tsv 2> {W}/filter_{g}.log",
                           shell=True, executable="/bin/bash") for g in GEN]
    rc_ = [p.wait() for p in ps]
    print("mapping exit codes", rc_, flush=True)
    if any(rc_):
        sys.exit(1)
    sh(f"cut -f1 {W}/pool_dep1.tsv {W}/pool_dep2.tsv | sort -u > {W}/pool.names")
    sh(f"/usr/local/bin/seqkit grep -j 8 -f {W}/pool.names {READS} -o {W}/pool.fq.gz")
    names = set(open(f"{W}/pool.names").read().split())
    by = collections.defaultdict(set)
    for g in GEN:
        for l in open(f"{W}/pool_{g}.tsv"):
            a, b = l.split(); by[b].add(a)
    v = sorted(len(by[t["id"]]) for t in tests)
    print("pool reads", len(names), "reads per fill: min", v[0], "median", v[len(v) // 2], "max", v[-1], flush=True)
    fished = set(x[1:].split()[0] for x in subprocess.run(f"zcat {FISHED} | awk 'NR%4==1'", shell=True, capture_output=True, text=True).stdout.splitlines())
    print("fished reads", len(fished), "of them in pool", len(fished & names), flush=True)
    k = by["dep1_gap0441"]
    print("TMED4 reads in pool", len(k), "contains ERR16751478.985790:", "ERR16751478.985790" in k, flush=True)
    print(sh(f"cat {W}/filter_dep1.log {W}/filter_dep2.log; ls -la {W}/pool.fq.gz"), flush=True)
    print("finished", flush=True)


def all_fills():
    global TH
    import multiprocessing as mp, pickle
    TH = 2
    os.makedirs(f"{W}/run", exist_ok=True)
    tests = pickle.load(open(f"{W}/tests.pkl", "rb"))
    by = collections.defaultdict(set)
    for g in GEN:
        for l in open(f"{W}/pool_{g}.tsv"):
            a, b = l.split(); by[b].add(a)
    p = subprocess.Popen(f"zcat {W}/pool.fq.gz", shell=True, stdout=subprocess.PIPE, text=True)
    while True:
        h = p.stdout.readline()
        if not h:
            break
        POOL[h.split()[0][1:]] = h + "".join(p.stdout.readline() for _ in range(3))
    print("pool in memory", len(POOL), flush=True)
    bad = [t for t in tests if not t["O"]]
    tests = [t for t in tests if t["O"]]
    for t in tests:
        t["names"] = sorted(by[t["id"]])
    out = open(f"{W}/fills.tsv", "w"); rr = open(f"{W}/fills.reads.tsv", "w")
    print(*HEAD, sep="\t", file=out)
    for t in bad:
        print(*t["meta"], "no_other_window", sep="\t", file=out)
    # guard: TMED4 must come out as the wrong allele, else stop before anything else runs
    k = next(t for t in tests if t["id"] == "dep1_gap0441")
    row, reads = one(k)
    print("TMED4", row, flush=True)
    if row[-1] != "wrong_allele":
        print("GUARD FAILED, stopping", flush=True)
        sys.exit(1)
    print(*row, sep="\t", file=out)
    for r in reads:
        print(*r, sep="\t", file=rr)
    n = 1
    with mp.Pool(12) as pl:
        for row, reads in pl.imap_unordered(one, [t for t in tests if t is not k]):
            print(*row, sep="\t", file=out); out.flush()
            for r in reads:
                print(*r, sep="\t", file=rr)
            n += 1
            if n % 25 == 0:
                print(n, "done", flush=True)
    print("finished", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "filter":
        filt(sys.argv[2])
        sys.exit(0)
    os.makedirs(W, exist_ok=True)
    if sys.argv[1] == "controls":
        controls()
    elif sys.argv[1] == "pool":
        pool()
    elif sys.argv[1] == "all":
        all_fills()
