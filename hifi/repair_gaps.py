#!/usr/bin/env python3
"""Repair the N gaps of a deposited haplotype with fills from the other haplotype / the asm25 haplotypes, verify every
new junction with HiFi reads, and write a main repaired assembly plus a 'gene-search' assembly.

Cause per gap (gap_classify.py label) -> action:
  hom_one_hap                  insert (A: two-copy stretch written into one haplotype; fill from the other deposited
                               haplotype, else an asm25 fill whose fill_cn call is two_copies)
  fill_nocall                  no fill from the other deposited haplotype: decided by the asm25 fill's fill_cn call,
                               two_copies -> insert (A); one_copy / unclear / short -> genesearch (B); repeat -> leave
  join                         join (C: flanks adjacent in another assembly; the <= 200 bp between them are inserted,
                               overlapping flanks are kept once)
  one_copy_fill, fill_unclear, fill_short_or_unscored
                               genesearch (B: insert in lower case into the gene-search assembly only)
  anything else                leave (D)
Splicing: the fill (fills.fa entry: fill + up to 1 kb of source contig on each side, in the gapped assembly's orientation)
is anchored on each side by the unique 25-mer of the deposited flank closest to the gap (unique within the 1.5 kb flank
window and within the fill entry, 50 bp exact match, position within 600 bp of the gapfill.py expectation). Deposited
bases between the anchor and the N run are replaced (counted as flank_bp_replaced).

Steps (OUT = per-assembly work dir):
  plan    repair_gaps.py plan GAPPED.fa GAPS.tsv FILLS.fa ALL_FILLS.fa FILL_CN.tsv CAUSES.tsv OTHER_DEP OUT [--gaps g,g]
          -> OUT/plan.tsv, OUT/pieces.fa (inserted sequences), OUT/windows.fa (+-15 kb around every repair, all repairs
             applied), OUT/junctions.tsv, OUT/cand.fa (25-mers from the 3 kb on both sides of every junction)
  words   repair_gaps.py words OUT COUNTS.txt     (jellyfish dump -c of cand.fa counted over GAPPED.fa + pieces.fa)
          -> OUT/anchors.txt (count == 1)
  verify  repair_gaps.py verify GAPPED.fa OUT READS.fq LABEL [THREADS]
          -> OUT/LABEL.repaired.fa, OUT/LABEL.genesearch.fa, OUT/repairs.tsv, OUT/summary.txt
          (--edited-only: write only chromosomes that carry an accepted repair; for toys)
  bed     repair_gaps.py bed OUT LABEL -> OUT/LABEL.repaired.bed, OUT/LABEL.genesearch.bed (inserted bases; a
          zero-length join is a 0-length interval)
A read crosses a junction when one alignment record (minimap2 -c map-hifi to windows.fa) covers >= 1 kb on both sides,
|read span - target span| <= max(200, 2 % of target span) and identity >= 0.97. A repair is accepted with >= 3 crossing
reads at every junction; otherwise it stays N (main) / N (gene-search) and the reason is recorded."""
import os
import re
import statistics
import subprocess
import sys
from collections import Counter, defaultdict

K, W, FO, SLACK = 25, 1500, 1000, 600
WIN, NEAR, STEP = 15000, 3000, 20
MIN_SIDE, MIN_READS, MIN_IDENT = 1000, 3, 0.97
comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
rc = lambda s: s.translate(comp)[::-1]
A_CAUSES = {"hom_one_hap"}
B_CAUSES = {"one_copy_fill", "fill_unclear", "fill_short_or_unscored"}


def read_fa(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for l in fh:
            if l.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = l[1:].split()[0], []
            else:
                buf.append(l.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def write_fa(path, seqs, width=80):
    with open(path, "w") as fh:
        for n, s in seqs.items():
            fh.write(f">{n}\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")


def anchor(g, a, b, lo_lim, hi_lim, P, fo_left, fill):
    """Splice coordinates: (x, y, eP, sP) -> new = g[:x] + P[eP:sP] + g[y:], or (None, reason)."""
    Pu = P.upper()
    pos = defaultdict(list)
    for i in range(len(Pu) - K + 1):
        pos[Pu[i:i + K]].append(i)
    fo_right = len(P) - fo_left - fill
    DL = g[max(lo_lim, a - W):a].upper()
    dl0 = a - len(DL)
    cL = Counter(DL[i:i + K] for i in range(len(DL) - K + 1))
    left = None
    for i in range(len(DL) - K, -1, -1):
        k = DL[i:i + K]
        if "N" in k or cL[k] != 1 or len(pos.get(k, ())) != 1:
            continue
        p = pos[k][0]
        if abs(p - (fo_left - (len(DL) - i))) > SLACK or i < 25 or p < 25 or DL[i - 25:i + K] != Pu[p - 25:p + K]:
            continue
        left = (dl0 + i + K, p + K)
        break
    if left is None:
        return None, "no_left_anchor"
    DR = g[b:min(hi_lim, b + W)].upper()
    cR = Counter(DR[j:j + K] for j in range(len(DR) - K + 1))
    for j in range(len(DR) - K + 1):
        k = DR[j:j + K]
        if "N" in k or cR[k] != 1 or len(pos.get(k, ())) != 1:
            continue
        p = pos[k][0]
        if p < left[1] or abs(p - (len(P) - fo_right + j)) > SLACK or DR[j:j + K + 25] != Pu[p:p + K + 25] \
                or len(DR[j:j + K + 25]) < K + 25:
            continue
        return (left[0], b + j, left[1], p), ""
    return None, "no_right_anchor"


def plan(argv):
    gfa, gtsv, ffa, allfa, cntsv, causes, other, out = argv[:8]
    only = set(argv[argv.index("--gaps") + 1].split(",")) if "--gaps" in argv else None
    os.makedirs(out, exist_ok=True)
    g = read_fa(gfa)
    rows = [l.rstrip("\n").split("\t") for l in open(gtsv)]
    h, rows = rows[0], rows[1:]
    ix = {k: i for i, k in enumerate(h)}
    labels = [c[:-7] for c in h if c.endswith("_status")]
    # fill_cn.tsv rows are in all_fills.fa order; keep the rows whose fill entry belongs to this assembly's fills.fa
    allh = [l[1:].strip() for l in open(allfa) if l.startswith(">")]
    cnr = [l.rstrip("\n").split("\t") for l in open(cntsv)][1:]
    assert len(allh) == len(cnr)
    fills = {}
    name = None
    for l in open(ffa):
        if l.startswith(">"):
            name = l[1:].strip(); fills[name] = []
        else:
            fills[name].append(l.strip())
    fills = {n: "".join(v) for n, v in fills.items()}
    call, fseq = {}, {}
    for hd, r in zip(allh, cnr):
        assert hd.split("|")[:2] == r[:2]
        if hd in fills:
            call[(r[0], r[1])] = r[6]
    for n, s in fills.items():
        f = n.split("|")
        lo = int(f[2].rsplit(":", 1)[1].split("-")[0])
        fill = int(f[4].split("=")[1])
        fol = min(FO, lo) if f[3] == "+" else len(s) - min(FO, lo) - fill   # '-': entry is reverse-complemented
        fseq[(f[0], f[1])] = (s, fol, fill, n)
    cause = {l.split("\t")[0]: l.rstrip("\n").split("\t")[-1] for l in list(open(causes))[1:]}
    runs = defaultdict(list)
    for r in rows:
        runs[r[ix["chrom"]]].append((int(r[ix["start"]]), int(r[ix["end"]])))
    for c in runs:
        runs[c].sort()
    order = [other] + [l for l in labels if l != other]

    def neighbours(c, a, b):
        rr = runs[c]
        k = rr.index((a, b))
        return (rr[k - 1][1] if k else 0), (rr[k + 1][0] if k + 1 < len(rr) else len(g[c]))

    plans = []
    for r in rows:
        gid, c, a, b = r[0], r[1], int(r[2]), int(r[3])
        if only and gid not in only:
            continue
        cs = cause.get(gid, "other")
        sp = [l for l in order if r[ix[f"{l}_status"]] == "spanned" and (gid, l) in fseq]
        if cs in A_CAUSES:
            act, cand = "insert", [l for l in sp if call.get((gid, l)) == "two_copies"]
        elif cs == "fill_nocall":
            cl = [call.get((gid, l)) for l in sp]
            if "two_copies" in cl:
                act, cand = "insert", [l for l in sp if call.get((gid, l)) == "two_copies"]
            elif cl and all(x == "repeat" for x in cl):
                act, cand = "leave", []
            else:
                act, cand = "genesearch", [l for l in sp if call.get((gid, l)) != "repeat"]
        elif cs == "join":
            act, cand = "join", [l for l in sp if fseq[(gid, l)][2] <= 200]
        elif cs in B_CAUSES:
            act, cand = "genesearch", sp
        else:
            act, cand = "leave", []
        rec = dict(gap=gid, chrom=c, start=a, end=b, cause=cs, action=act, source="", src_call="", fill=0,
                   x=a, y=b, ins="", reason="", repl_L=0, repl_R=0)
        if act != "leave":
            why = "no_source"
            lim = neighbours(c, a, b)
            for l in cand:
                P, fol, fill, hd = fseq[(gid, l)]
                if re.search("[Nn]", P):
                    why = "fill_has_N"; continue
                sp_, why = anchor(g[c], a, b, lim[0], lim[1], P, fol, fill)
                if sp_ is None:
                    continue
                x, y, eP, sP = sp_
                ins = P[eP:sP].upper()
                rec.update(source=hd, src_call=call.get((gid, l), ""), fill=fill, x=x, y=y,
                           ins=ins.lower() if act == "genesearch" else ins, repl_L=a - x, repl_R=y - b)
                why = ""
                break
            if why:
                rec["reason"] = why
                rec["action"] = act + ":not_built"
        plans.append(rec)
    # candidate sequence with every built repair; windows and junctions in its coordinates
    by = defaultdict(list)
    for p in plans:
        if p["action"] in ("insert", "join", "genesearch"):
            by[p["chrom"]].append(p)
    wins, junc, cand = {}, [], set()
    for c, pp in by.items():
        pp.sort(key=lambda p: p["x"])
        s, parts, prev, shift = g[c], [], 0, 0
        for p in pp:
            assert p["x"] >= prev, (p["gap"], "overlapping edits")
            parts.append(s[prev:p["x"]])
            p["nx"] = p["x"] + shift
            p["ny"] = p["nx"] + len(p["ins"])
            parts.append(p["ins"])
            shift += len(p["ins"]) - (p["y"] - p["x"])
            prev = p["y"]
        parts.append(s[prev:])
        ns = "".join(parts)
        iv = []
        for p in pp:
            iv.append([max(0, p["nx"] - WIN), min(len(ns), p["ny"] + WIN), [p]])
        merged = []
        for v in iv:
            if merged and v[0] <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], v[1]); merged[-1][2] += v[2]
            else:
                merged.append(v)
        for ws, we, pl in merged:
            wn = f"{c}:{ws}-{we}"
            wins[wn] = ns[ws:we]
            for p in pl:
                js = [("L", p["nx"]), ("R", p["ny"])] if p["ny"] > p["nx"] else [("J", p["nx"])]
                for side, j in js:
                    junc.append((p["gap"], side, wn, j - ws))
                    for t in (ns[max(0, j - NEAR):j], ns[j:j + NEAR]):
                        t = t.upper()
                        for i in range(0, len(t) - K + 1, STEP):
                            k = t[i:i + K]
                            if "N" not in k:
                                cand.add(min(k, rc(k)))
    write_fa(os.path.join(out, "windows.fa"), wins)
    with open(os.path.join(out, "pieces.fa"), "w") as fh:
        for p in plans:
            if p["ins"]:
                fh.write(f">{p['gap']}\n{p['ins'].upper()}\n")
    with open(os.path.join(out, "cand.fa"), "w") as fh:
        fh.write("".join(f">{i}\n{k}\n" for i, k in enumerate(sorted(cand))))
    with open(os.path.join(out, "junctions.tsv"), "w") as fh:
        fh.write("gap\tside\twindow\tpos\n" + "".join("\t".join(map(str, j)) + "\n" for j in junc))
    cols = ["gap", "chrom", "start", "end", "cause", "action", "source", "src_call", "fill", "x", "y", "repl_L",
            "repl_R", "ins_len", "reason"]
    with open(os.path.join(out, "plan.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for p in plans:
            p["ins_len"] = len(p["ins"])
            fh.write("\t".join(str(p[k]) for k in cols) + "\n")
    print(Counter(p["action"] for p in plans), f"{len(wins)} windows, {len(junc)} junctions, {len(cand)} words",
          file=sys.stderr)


def words(argv):
    out, counts = argv[:2]
    keep = [k for k, n in (x.split() for x in open(counts)) if n == "1"]
    open(os.path.join(out, "anchors.txt"), "w").write("\n".join(keep) + "\n")
    print(f"{len(keep)} unique words", file=sys.stderr)


def verify(argv):
    gfa, out, reads, label = argv[:4]
    threads = argv[4] if len(argv) > 4 and not argv[4].startswith("--") else "16"
    edited_only = "--edited-only" in argv          # toy runs: write only the chromosomes that carry a repair
    junc = [l.rstrip("\n").split("\t") for l in open(os.path.join(out, "junctions.tsv"))][1:]
    jby = defaultdict(list)
    for gid, side, wn, p in junc:
        jby[wn].append((gid, side, int(p)))
    cross = Counter()
    seen = set()
    pr = subprocess.Popen(["minimap2", "-c", "-x", "map-hifi", "--secondary=no", "-t", threads,
                           os.path.join(out, "windows.fa"), reads], stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL, text=True)
    for l in pr.stdout:
        f = l.split("\t")
        q, qs, qe, wn, ts, te, nm, bl = f[0], int(f[2]), int(f[3]), f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])
        if nm / bl < MIN_IDENT or abs((qe - qs) - (te - ts)) > max(200, 0.02 * (te - ts)):
            continue
        for gid, side, j in jby.get(wn, ()):
            if ts <= j - MIN_SIDE and te >= j + MIN_SIDE and (q, gid, side) not in seen:
                seen.add((q, gid, side))
                cross[(gid, side)] += 1
    if pr.wait():
        sys.exit("minimap2 failed")
    pl = [l.rstrip("\n").split("\t") for l in open(os.path.join(out, "plan.tsv"))]
    ph, pl = pl[0], pl[1:]
    ix = {k: i for i, k in enumerate(ph)}
    pieces = read_fa(os.path.join(out, "pieces.fa"))
    sides = defaultdict(list)
    for gid, side, _, _ in junc:
        sides[gid].append(side)
    g = read_fa(gfa)
    rep = []
    for r in pl:
        gid, act = r[ix["gap"]], r[ix["action"]]
        jl = [cross[(gid, s)] for s in sides.get(gid, [])]
        if act in ("insert", "join", "genesearch"):
            ok = len(jl) > 0 and min(jl) >= MIN_READS
            why = "" if ok else f"junction reads < {MIN_READS}"
        else:
            ok, why = False, r[ix["reason"]] or ("not repairable (leave N)" if act == "leave" else "")
        jL = jl[0] if jl else ""
        jR = jl[-1] if jl else ""
        rep.append(dict(gap=gid, chrom=r[ix["chrom"]], start=int(r[ix["start"]]), end=int(r[ix["end"]]),
                        cause=r[ix["cause"]], action=act, ins_len=int(r[ix["ins_len"]]), source=r[ix["source"]],
                        src_call=r[ix["src_call"]], jL=jL, jR=jR, ok=ok, reason=why, x=int(r[ix["x"]]),
                        y=int(r[ix["y"]]), repl=f"{r[ix['repl_L']]}/{r[ix['repl_R']]}"))
    by_main, by_gs = defaultdict(list), defaultdict(list)
    for p in rep:
        if p["ok"]:
            seq = pieces.get(p["gap"], "")          # zero-length joins have no pieces.fa entry
            if p["action"] == "genesearch":
                by_gs[p["chrom"]].append((p["x"], p["y"], seq.lower()))
            else:
                by_main[p["chrom"]].append((p["x"], p["y"], seq))
                by_gs[p["chrom"]].append((p["x"], p["y"], seq))

    if edited_only:
        g = {c: s for c, s in g.items() if c in by_gs or c in by_main}

    def apply(edits):
        res = {}
        for c, s in g.items():
            ee = sorted(edits.get(c, []))
            if not ee:
                res[c] = s; continue
            parts, prev = [], 0
            for x, y, ins in ee:
                parts += [s[prev:x], ins]; prev = y
            parts.append(s[prev:])
            res[c] = "".join(parts)
        return res

    nrun = lambda seqs: sum(len(re.findall("[Nn]{10,}", s)) for s in seqs.values())
    main = apply(by_main)
    write_fa(os.path.join(out, f"{label}.repaired.fa"), main)
    n_main, len_main = nrun(main), sum(map(len, main.values()))
    del main
    gs = apply(by_gs)
    write_fa(os.path.join(out, f"{label}.genesearch.fa"), gs)
    n_gs, len_gs = nrun(gs), sum(map(len, gs.values()))
    del gs
    cols = ["assembly", "gap", "chrom", "start", "end", "cause", "action", "ins_len", "flank_bp_replaced_LR", "source",
            "src_call", "junction_reads_L", "junction_reads_R", "accepted", "reason"]
    with open(os.path.join(out, "repairs.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for p in rep:
            fh.write("\t".join(map(str, [label, p["gap"], p["chrom"], p["start"], p["end"], p["cause"], p["action"],
                                         p["ins_len"], p["repl"], p["source"], p["src_call"], p["jL"], p["jR"],
                                         "yes" if p["ok"] else "no", p["reason"]])) + "\n")
    L = []
    L.append(f"{label}: N runs >= 10 before {nrun(g)}, after: main {n_main}, gene-search {n_gs}")
    L.append(f"{label}: length before {sum(map(len, g.values()))}, main {len_main}, gene-search {len_gs}")
    L.append("action\tcause\tn\taccepted\trejected\tbp_inserted_accepted\tmedian_min_junction_reads")
    grp = defaultdict(list)
    for p in rep:
        grp[(p["action"], p["cause"])].append(p)
    for (act, cs), pp in sorted(grp.items()):
        acc = [p for p in pp if p["ok"]]
        mj = [min(p["jL"], p["jR"]) for p in pp if p["jL"] != ""]
        L.append(f"{act}\t{cs}\t{len(pp)}\t{len(acc)}\t{len(pp) - len(acc)}\t{sum(p['ins_len'] for p in acc)}\t"
                 f"{statistics.median(mj) if mj else '-'}")
    for act in ("insert", "join", "genesearch"):
        pp = [p for p in rep if p["action"] == act]
        mj = [min(p["jL"], p["jR"]) for p in pp]
        L.append(f"{act}: {len(pp)} built, {sum(p['ok'] for p in pp)} accepted, "
                 f"{sum(p['ins_len'] for p in pp if p['ok'])} bp inserted, median crossing reads (lower junction) "
                 f"{statistics.median(mj) if mj else '-'}")
    open(os.path.join(out, "summary.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


def bed(argv):
    """Coordinates of the accepted repairs in the output FASTAs (the deposited assemblies are soft-masked, so the lower
    case of gene-search inserts alone does not mark them): OUT/LABEL.repaired.bed, OUT/LABEL.genesearch.bed."""
    out, label = argv[:2]
    pl = [l.rstrip("\n").split("\t") for l in open(os.path.join(out, "plan.tsv"))]
    ix = {k: i for i, k in enumerate(pl[0])}
    xy = {r[0]: (r[ix["chrom"]], int(r[ix["x"]]), int(r[ix["y"]]), int(r[ix["ins_len"]])) for r in pl[1:]}
    rp = [l.rstrip("\n").split("\t") for l in open(os.path.join(out, "repairs.tsv"))]
    jx = {k: i for i, k in enumerate(rp[0])}
    acc = [r for r in rp[1:] if r[jx["accepted"]] == "yes"]
    for kind in ("repaired", "genesearch"):
        ed = defaultdict(list)
        for r in acc:
            if kind == "repaired" and r[jx["action"]] == "genesearch":
                continue
            c, x, y, n = xy[r[jx["gap"]]]
            ed[c].append((x, y, n, f"{r[jx['gap']]}|{r[jx['cause']]}|{r[jx['action']]}"))
        with open(os.path.join(out, f"{label}.{kind}.bed"), "w") as fh:
            for c in sorted(ed):
                shift = 0
                for x, y, n, name in sorted(ed[c]):
                    fh.write(f"{c}\t{x + shift}\t{x + shift + n}\t{name}\n")
                    shift += n - (y - x)


if __name__ == "__main__":
    {"plan": plan, "words": words, "verify": verify, "bed": bed}[sys.argv[1]](sys.argv[2:])
