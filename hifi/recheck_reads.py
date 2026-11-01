#!/usr/bin/env python3
"""Read-level re-examination of gap labels on the deposited hap1 windows (recheck_prep.py win.fa).

Inputs: WIN.fa (h1 windows `gap|h1|chrom:a-b|gap=x-y` and alt windows `gap|alt|src|contig:a-b|fill=x-y`), H1.paf and
ALT.paf (minimap2 -c --cs -x map-hifi of the fished reads to the h1 windows / to the alt windows, separately), READS.fq
(for tails), WIN.trf.dat (trf -ngs on WIN.fa), OUT dir.
Per h1 window:
  het sites: >= 6 reads, alt base in 25-75 % (alignments >= 2 kb). A read is h1 (carries the deposited hap1 base at
            >= 80 % of >= 2 informative sites), h2 (alt base at >= 80 %) or ? (else / too few sites).
  gap-end events (side L: alignment reaching within EDGE of the N run from the left, anchored >= 1 kb; side R mirror):
            span      the same read aligns on the other side of the N run, same strand, read order consistent;
                      size = read bases between the two N-run edges (the real content of the gap on that molecule)
            fill      the read continues >= 1 kb and its alt-window alignment covers the alt fill edge (copy carrying
                      the fill)
            clip      continues >= 1 kb with sequence aligned to neither (tail written to tails.fa)
            end       fewer than 1 kb of read left
  approach  alignments that stop 0.3-15 kb before the N run (toward it): clip (read goes on >= 1 kb) / end
  depth     reads per 500 bp bin (primary alignments >= 2 kb), total / h1 / h2
  tandem    TRF arrays (>= 200 bp, merged within 300 bp) touching the gap (within 2 kb) on each side: period, % match,
            array length in the assembly; reads that enter it from unique sequence (alignment starts >= 500 bp before
            the array) followed through alignment + tail: span (reaches the far side), exit (tail leaves array-like
            sequence), in_array_end (read ends still in array-like sequence); array bases seen per read
  tract     GA/TC tract (>= 100 bp, >= 90 % purine or pyrimidine) touching the gap: depth in tract vs 2-10 kb away, reads
            through the tract in one alignment, their error rate inside vs outside the tract, alignments ending
            inside the tract
Writes OUT/events.tsv, OUT/summary.tsv, OUT/depth.tsv, OUT/tandem.tsv, OUT/tract.tsv, OUT/tails.fa, OUT/readlen.txt
Usage: recheck_reads.py WIN.fa H1.paf ALT.paf READS.fq WIN.trf.dat OUT"""
import os
import re
import statistics
import sys
from collections import Counter, defaultdict

EDGE, ANCH, MIN_CONT, TAIL = 300, 1000, 1000, 6000
APPROACH = 15000
comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
rc = lambda s: s.translate(comp)[::-1]
tok = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+)")
winfa, h1paf, altpaf, fq, trf, out = sys.argv[1:7]
os.makedirs(out, exist_ok=True)

seqs, name = {}, None
for l in open(winfa):
    if l.startswith(">"):
        name = l[1:].strip(); seqs[name] = []
    else:
        seqs[name].append(l.strip().upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
h1 = {}
alt = {}
for n in seqs:
    f = n.split("|")
    if f[1] == "h1":
        x, y = map(int, f[3].split("=")[1].split("-"))
        h1[f[0]] = (n, x, y)
    else:
        x, y = map(int, f[4].split("=")[1].split("-"))
        alt[f[0]] = (n, x, y)


def load(paf):
    al = defaultdict(list)          # target -> list of records
    byread = defaultdict(list)
    for l in open(paf):
        f = l.rstrip("\n").split("\t")
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), "")
        tp = next((x[5:] for x in f[12:] if x.startswith("tp:A:")), "P")
        r = dict(q=f[0], ql=int(f[1]), qs=int(f[2]), qe=int(f[3]), st=f[4], t=f[5], ts=int(f[7]), te=int(f[8]),
                 nm=int(f[9]), bl=int(f[10]), cs=cs, tp=tp)
        al[r["t"]].append(r)
        byread[(r["q"], r["t"])].append(r)
    return al, byread


A1, B1 = load(h1paf)
A2, B2 = load(altpaf)
rlen = {}
need_tail = set()


def walk(r):
    """yield (kind, target_pos, length, base) along the cs string"""
    p = r["ts"]
    for m in tok.findall(r["cs"]):
        if m[0] == ":":
            n = int(m[1:]); yield ("=", p, n, None); p += n
        elif m[0] == "*":
            yield ("x", p, 1, m[2].upper()); p += 1
        elif m[0] == "+":
            yield ("i", p, len(m) - 1, None)
        else:
            yield ("d", p, len(m) - 1, None); p += len(m) - 1


ev_out = open(os.path.join(out, "events.tsv"), "w")
ev_out.write("gap\tside\tread\thap\tstrand\tts\tte\tqlen\tevent\tvalue\n")
sm_out = open(os.path.join(out, "summary.tsv"), "w")
sm_out.write("gap\tn_het\tside\thap\treach\tspan\tfill\tclip\tend\tspan_sizes\tapproach_clip\tapproach_end\n")
dp_out = open(os.path.join(out, "depth.tsv"), "w")
dp_out.write("gap\tbin_start\tdist_to_gap\tall\th1\th2\n")
tn_out = open(os.path.join(out, "tandem.tsv"), "w")
tn_out.write("gap\tside\tarr_start\tarr_end\tarr_len\tperiod\tpct_match\tcopies\tn_enter\tspan\texit\tin_array_end\t"
             "array_bp_in_reads(max,median)\tread_rest_after_entry(median)\tspan_array_len\n")
tr_out = open(os.path.join(out, "tract.tsv"), "w")
tr_out.write("gap\tside\ttract_start\ttract_end\ttract_len\tpurity\tdepth_tract\tdepth_flank\tthrough\tthrough_err_in\t"
             "through_err_out\tend_in_tract\tclip_in_tract\n")
tails = {}

# TRF arrays per window
arrays = defaultdict(list)
cur = None
for l in open(trf):
    if l.startswith("@"):
        cur = l[1:].strip()
    elif cur and l[:1].isdigit():
        f = l.split()
        a, b = int(f[0]) - 1, int(f[1])
        if b - a >= 200:
            arrays[cur].append((a, b, int(f[2]), float(f[3]), int(f[5]), f[13] if len(f) > 13 else ""))


def merged(arrs):
    arrs = sorted(arrs)
    m = []
    for a in arrs:
        if m and a[0] <= m[-1][1] + 300:
            m[-1][1] = max(m[-1][1], a[1]); m[-1][2].append(a)
        else:
            m.append([a[0], a[1], [a]])
    return m


def kset(s, k=15):
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def frac_in(s, ks, k=15):
    n = len(s) - k + 1
    if n <= 0:
        return 0
    return sum((s[i:i + k] in ks) or (rc(s[i:i + k]) in ks) for i in range(0, n, 3)) / len(range(0, n, 3))


all_lens = []
results = []
for g, (wn, x, y) in sorted(h1.items()):
    recs = [r for r in A1.get(wn, []) if r["tp"] in "PS"]
    W = len(seqs[wn])
    # het sites
    cov = defaultdict(int)
    subs = defaultdict(Counter)
    for r in recs:
        if r["te"] - r["ts"] < 2000:
            continue
        for b in range(r["ts"] // 500, r["te"] // 500 + 1):
            cov[b] += 1
        for kind, p, n, base in walk(r):
            if kind == "x":
                subs[p][base] += 1
    het = {}
    for p, d in subs.items():
        base, n = d.most_common(1)[0]
        c = cov[p // 500]
        if c >= 6 and n >= 2 and 0.25 <= n / c <= 0.75:
            het[p] = base
    hap = {}
    reads = defaultdict(list)
    for r in recs:
        reads[r["q"]].append(r)
    for q, rr in reads.items():
        ref = alt_ = 0
        for r in rr:
            xs = {p: b for k, p, n, b in walk(r) if k == "x"}
            for p in het:
                if r["ts"] <= p < r["te"]:
                    if xs.get(p) == het[p]:
                        alt_ += 1
                    elif p not in xs:
                        ref += 1
        tot = ref + alt_
        hap[q] = "h1" if tot >= 2 and ref / tot >= 0.8 else ("h2" if tot >= 2 and alt_ / tot >= 0.8 else "?")
    # depth
    bins = defaultdict(Counter)
    for r in recs:
        if r["tp"] != "P" or r["te"] - r["ts"] < 2000:
            continue
        for b in range(r["ts"] // 500, (r["te"] - 1) // 500 + 1):
            bins[b]["all"] += 1; bins[b][hap[r["q"]]] += 1
    for b in range(0, W // 500 + 1):
        s0 = b * 500
        d = (s0 - x) if s0 < x else (s0 - y if s0 >= y else 0)
        dp_out.write(f"{g}\t{s0}\t{d}\t{bins[b]['all']}\t{bins[b]['h1']}\t{bins[b]['h2']}\n")
    # gap-end events
    an = alt.get(g)
    cnt = defaultdict(Counter)
    sizes = defaultdict(list)
    seen = set()
    for q, rr in reads.items():
        ql = rr[0]["ql"]
        for r in rr:
            if r["ts"] <= x - ANCH and r["te"] >= y + ANCH:      # one alignment straight through the N run
                for side in "LR":
                    if (q, side) not in seen:
                        seen.add((q, side))
                        val = (y - x) + ((r["qe"] - r["qs"]) - (r["te"] - r["ts"]))
                        cnt[(side, hap[q])]["span"] += 1
                        sizes[(side, hap[q])].append(val)
                        ev_out.write(f"{g}\t{side}\t{q}\t{hap[q]}\t{r['st']}\t{r['ts']}\t{r['te']}\t{ql}\tspan\t{val}\n")
                continue
            for side in "LR":
                if side == "L":
                    reach = r["te"] >= x - EDGE and r["ts"] <= x - ANCH and r["te"] <= y
                    near = (not reach) and x - APPROACH <= r["te"] < x - EDGE and r["ts"] < r["te"] - ANCH
                    fwd = r["st"] == "+"          # moving rightwards in the window = forward in the read
                else:
                    reach = r["ts"] <= y + EDGE and r["te"] >= y + ANCH and r["ts"] >= x
                    near = (not reach) and y + EDGE < r["ts"] <= y + APPROACH and r["te"] > r["ts"] + ANCH
                    fwd = r["st"] == "-"          # moving leftwards in the window = forward in the read for '-'
                rest = (ql - r["qe"]) if fwd else r["qs"]
                if near:
                    ev = "approach_clip" if rest >= MIN_CONT else "approach_end"
                    cnt[(side, hap[q])][ev] += 1
                    ev_out.write(f"{g}\t{side}\t{q}\t{hap[q]}\t{r['st']}\t{r['ts']}\t{r['te']}\t{ql}\t{ev}\t{rest}\n")
                    continue
                if not reach or (q, side) in seen:
                    continue
                seen.add((q, side))
                ev, val = None, ""
                for o in rr:
                    if o is r or o["st"] != r["st"]:
                        continue
                    if side == "L" and o["ts"] >= y - EDGE and o["te"] > y + 500:
                        if r["st"] == "+" and o["qs"] >= r["qe"] - 500:
                            ev, val = "span", (o["qs"] - (o["ts"] - y)) - (r["qe"] + (x - r["te"]))
                        elif r["st"] == "-" and o["qe"] <= r["qs"] + 500:
                            ev, val = "span", (r["qs"] - (x - r["te"])) - (o["qe"] + (o["ts"] - y))
                    if side == "R" and o["te"] <= x + EDGE and o["ts"] < x - 500:
                        if r["st"] == "+" and o["qe"] <= r["qs"] + 500:
                            ev, val = "span", (r["qs"] - (r["ts"] - y)) - (o["qe"] + (x - o["te"]))
                        elif r["st"] == "-" and o["qs"] >= r["qe"] - 500:
                            ev, val = "span", (o["qs"] - (x - o["te"])) - (r["qe"] + (r["ts"] - y))
                if ev is None:
                    if rest < MIN_CONT:
                        ev, val = "end", rest
                    else:
                        ev, val = "clip", rest
                        if an:
                            an_n, fx, fy = an
                            edge = fx if side == "L" else fy
                            for o in B2.get((q, an_n), []):
                                if o["ts"] <= edge - 500 and o["te"] >= edge + 500:
                                    ev = "fill"; break
                        if ev == "clip":
                            need_tail.add(q)
                            tails[(g, side, q)] = (r["st"], r["qs"], r["qe"], fwd)
                cnt[(side, hap[q])][ev] += 1
                if ev == "span":
                    sizes[(side, hap[q])].append(val)
                ev_out.write(f"{g}\t{side}\t{q}\t{hap[q]}\t{r['st']}\t{r['ts']}\t{r['te']}\t{ql}\t{ev}\t{val}\n")
    for side in "LR":
        for hp in ("h1", "h2", "?"):
            c = cnt[(side, hp)]
            reach = sum(c[k] for k in ("span", "fill", "clip", "end"))
            sm_out.write(f"{g}\t{len(het)}\t{side}\t{hp}\t{reach}\t{c['span']}\t{c['fill']}\t{c['clip']}\t{c['end']}\t"
                         f"{','.join(map(str, sorted(sizes[(side, hp)])))}\t{c['approach_clip']}\t{c['approach_end']}\n")
    results.append((g, wn, x, y, reads, hap))

# reads for tails / tandem tails
need = set(need_tail)
for g, wn, x, y, reads, hap in results:
    need |= set(reads)
rseq = {}
with open(fq) as fh:
    while True:
        h = fh.readline()
        if not h:
            break
        s = fh.readline().strip(); fh.readline(); fh.readline()
        q = h[1:].split()[0]
        all_lens.append(len(s))
        if q in need:
            rseq[q] = s
with open(os.path.join(out, "readlen.txt"), "w") as fh:
    L = sorted(all_lens)
    fh.write(f"n {len(L)} mean {statistics.mean(L):.0f} median {L[len(L) // 2]} n50 "
             f"{next(v for v, c in zip(L[::-1], __import__('itertools').accumulate(L[::-1])) if c >= sum(L) / 2)}\n")
    fh.write(" ".join(map(str, L)) + "\n")
with open(os.path.join(out, "tails.fa"), "w") as fh:
    for (g, side, q), (st, qs, qe, fwd) in tails.items():
        s = rseq.get(q)
        if not s:
            continue
        t = s[qe:qe + TAIL] if fwd else rc(s[max(0, qs - TAIL):qs])
        fh.write(f">{g}|{side}|{q}\n{t}\n")


for g, wn, x, y, reads, hap in results:
    S = seqs[wn]
    m = merged(arrays.get(wn, []))
    # tandem arrays touching the gap
    for side in "LR":
        blk = None
        for a, b, parts in m:
            if side == "L" and b >= x - 2000 and a < x:
                blk = (a, min(b, x), parts)
            if side == "R" and a <= y + 2000 and b > y and blk is None:
                blk = (max(a, y), b, parts)
        if not blk:
            continue
        a, b, parts = blk
        top = max(parts, key=lambda p: p[1] - p[0])
        ks = kset(S[a:b])
        n_enter = span = ex = inend = 0
        arr_bp, rests, span_len = [], [], []
        for q, rr in reads.items():
            s = rseq.get(q)
            if not s:
                continue
            for r in rr:
                if side == "L":
                    ok = r["ts"] <= a - 500 and r["te"] >= a + 200
                    fwd = r["st"] == "+"
                else:
                    ok = r["te"] >= b + 500 and r["ts"] <= b - 200
                    fwd = r["st"] == "-"
                if not ok:
                    continue
                n_enter += 1
                # read position of the array entry point
                if side == "L":
                    pos_entry = (r["qs"] + (a - r["ts"])) if r["st"] == "+" else (r["qe"] - (a - r["ts"]))
                else:
                    pos_entry = (r["qe"] - (r["te"] - b)) if r["st"] == "+" else (r["qs"] + (r["te"] - b))
                rest = (len(s) - pos_entry) if fwd else pos_entry
                rests.append(rest)
                # does the read have an alignment on the far side of the gap?
                far = [o for o in rr if o is not r and o["st"] == r["st"] and
                       ((side == "L" and o["ts"] >= y - EDGE and o["te"] > y + 500) or
                        (side == "R" and o["te"] <= x + EDGE and o["ts"] < x - 500))]
                if far:
                    span += 1
                    o = far[0]
                    if side == "L":
                        p2 = (o["qs"] - (o["ts"] - y)) if r["st"] == "+" else (o["qe"] + (o["ts"] - y))
                    else:
                        p2 = (o["qe"] + (x - o["te"])) if r["st"] == "+" else (o["qs"] - (x - o["te"]))
                    span_len.append(abs(p2 - pos_entry))
                    break
                # follow array-like sequence from the entry point
                t = s[pos_entry:] if fwd else rc(s[:pos_entry])
                inarr, exited = 0, False
                for i in range(0, len(t), 500):
                    ch = t[i:i + 500]
                    if len(ch) < 100:
                        break
                    if frac_in(ch, ks) >= 0.3:
                        inarr = i + len(ch)
                    else:
                        exited = True
                        break
                arr_bp.append(inarr)
                if exited:
                    ex += 1
                else:
                    inend += 1
                break
        ab = f"{max(arr_bp)},{int(statistics.median(arr_bp))}" if arr_bp else ""
        tn_out.write(f"{g}\t{side}\t{a}\t{b}\t{b - a}\t{top[2]}\t{top[4]}\t{top[3]:.1f}\t{n_enter}\t{span}\t{ex}\t{inend}\t"
                     f"{ab}\t{int(statistics.median(rests)) if rests else ''}\t{','.join(map(str, sorted(span_len)))}\n")
    # GA/TC tract touching the gap
    for side in "LR":
        reg = (max(0, x - 3000), x) if side == "L" else (y, min(len(S), y + 3000))
        sub = S[reg[0]:reg[1]]
        best = None
        for alph in ("AG", "CT"):
            i = 0
            while i < len(sub):
                j, bad = i, 0
                while j < len(sub) and (bad + (sub[j] not in alph)) <= 0.1 * (j - i + 1):
                    bad += sub[j] not in alph; j += 1
                if j - i >= 100 and (best is None or j - i > best[1] - best[0]):
                    best = (i, j, alph, 1 - bad / max(1, j - i))
                i = max(j, i + 1)
        if not best:
            continue
        ta, tb = reg[0] + best[0], reg[0] + best[1]
        recs = [r for rr in reads.values() for r in rr if r["tp"] == "P"]
        dcount = lambda p: sum(1 for r in recs if r["ts"] <= p < r["te"])
        dt = statistics.mean(dcount(p) for p in range(ta, tb, 25))
        fl = (ta - 10000, ta - 2000) if side == "L" else (tb + 2000, tb + 10000)
        df = statistics.mean(dcount(p) for p in range(max(0, fl[0]), max(1, fl[1]), 500))
        thr, ein, eout, endin, clipin = 0, [], [], 0, 0
        for r in recs:
            if r["ts"] <= ta - 200 and r["te"] >= tb + 200:
                thr += 1
                e_in = e_out = 0
                for kind, p, n, base in walk(r):
                    if kind == "=":
                        continue
                    if ta <= p < tb:
                        e_in += n
                    elif ta - 2000 <= p < ta or tb <= p < tb + 2000:
                        e_out += n
                ein.append(1000 * e_in / (tb - ta)); eout.append(1000 * e_out / 4000)
            else:
                for endp, fwd in ((r["te"], r["st"] == "+"), (r["ts"], r["st"] == "-")):
                    if ta <= endp <= tb:
                        rest = (r["ql"] - r["qe"]) if fwd else r["qs"]
                        if rest >= MIN_CONT:
                            clipin += 1
                        else:
                            endin += 1
        md = lambda v: f"{statistics.median(v):.1f}" if v else ""
        tr_out.write(f"{g}\t{side}\t{ta}\t{tb}\t{tb - ta}\t{best[2]}:{best[3]:.2f}\t{dt:.1f}\t{df:.1f}\t{thr}\t{md(ein)}\t{md(eout)}\t{endin}\t{clipin}\n")
print("done", len(results), "windows", file=sys.stderr)
