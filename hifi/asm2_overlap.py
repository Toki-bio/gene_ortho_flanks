#!/usr/bin/env python3
"""Overlap merges: N runs whose two contig ends carry the SAME sequence (typically a tandem array written on both sides
of the N run; read-measured gap content is negative, recheck_label.py). The repair drops the N run and one copy of
the duplicated stretch: new = g[:a] + g[b + o:], o = overlap length. The result is written as a repair_gaps.py plan
directory (action 'join', ins_len 0) so that 'repair_gaps.py words' and 'verify' (>= 3 reads crossing the merged
junction with >= 1 kb each side, identity >= 0.97, span agreement max(200 bp, 2 %)) decide acceptance unchanged.

Overlap search: L = up to E + PAD bp before the N run, R = up to E + PAD bp after it (never across a neighbouring N
run). Every exact K-mer of L's last TAIL bp that occurs in R votes for o = len(L) - i + j (L[i] <-> R[j]); K-mers
occurring > MAXOCC times in R are skipped. Best o = most votes within +-SMOOTH bp; when E (expected overlap from reads)
is given, o is searched in E +- max(MINWIN, 0.3 E), otherwise in [MINO, len(R)]. ident = share of the overlap's
K-mer positions (whole overlap, +-SMOOTH diagonals) that match; o is kept when ident >= MIN_FRAC (a true duplicate is
~1.0; array-period shifts at 80-93 % unit identity give <= ~0.3).

If no exact duplicate is found and E plus the tandem-array lengths at the two ends (aL, aR; recheck tandem.tsv) are
given, the deposited ends are taken to carry E bp more array than the reads ("array over-length"): E bp of array are
removed, min(E, aR) from the right end's array next to the N run and the rest (<= aL) from the left end's array.

Usage: asm2_overlap.py GAPPED.fa GAPS.tsv OUT_DIR GAP[:E[:aL:aR]][,...]
  GAPS.tsv: gapfill.py gaps.tsv of GAPPED.fa (gap, chrom, start, end). E: expected overlap in bp (positive), optional.
  -> OUT_DIR/plan.tsv, pieces.fa (empty), windows.fa, junctions.tsv, cand.fa, overlap.tsv (search details)"""
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repair_gaps import read_fa, write_fa, rc, K as WK, WIN, NEAR, STEP  # noqa: E402

K, PAD, TAIL, MAXOCC, SMOOTH, MINWIN, MINO, MIN_FRAC, MAXO = 21, 6000, 3000, 50, 10, 3000, 200, 0.8, 30000

gfa, gtsv, out, spec = sys.argv[1:5]
os.makedirs(out, exist_ok=True)
want, arr = {}, {}
for x in spec.split(","):
    f = x.split(":")
    want[f[0]] = float(f[1]) if len(f) > 1 and f[1] else None
    if len(f) > 3:
        arr[f[0]] = (int(f[2]), int(f[3]))
g = read_fa(gfa)
rows = [l.rstrip("\n").split("\t") for l in open(gtsv)][1:]
runs = defaultdict(list)
for r in rows:
    runs[r[1]].append((int(r[2]), int(r[3])))
for c in runs:
    runs[c].sort()


def kmers(s):
    s = s.upper()
    return [(i, s[i:i + K]) for i in range(len(s) - K + 1) if "N" not in s[i:i + K]]


def search(L, R, E):
    rpos = defaultdict(list)
    for j, k in kmers(R):
        rpos[k].append(j)
    if E:
        lo, hi = max(MINO, E - max(MINWIN, 0.3 * E)), E + max(MINWIN, 0.3 * E)
    else:
        lo, hi = MINO, len(R)
    hi = min(hi, len(R), len(L))
    votes = Counter()
    for i, k in kmers(L[-TAIL:]):
        i += len(L) - min(TAIL, len(L))
        js = rpos.get(k, ())
        if len(js) > MAXOCC:
            continue
        for j in js:
            o = len(L) - i + j
            if lo <= o <= hi:
                votes[o] += 1
    if not votes:
        return None
    sm = {o: sum(votes.get(o + d, 0) for d in range(-SMOOTH, SMOOTH + 1)) for o in votes}
    best = max(sm, key=lambda o: (sm[o], votes[o]))
    o = max(range(best - SMOOTH, best + SMOOTH + 1), key=lambda x: votes.get(x, 0))
    # identity over the whole overlap: L[len(L)-o+t] vs R[t], t < o, exact K-mers on diagonals o +- SMOOTH
    Lk = {}
    for i, k in kmers(L[len(L) - o:]):
        Lk.setdefault(k, []).append(i)
    hit = set()
    for j, k in kmers(R[:o]):
        for i in Lk.get(k, ()):
            if abs(i - j) <= SMOOTH:
                hit.add(j); break
    frac = len(hit) / max(1, o - K + 1)
    return o, votes[o], sm[best], frac


plans, info = [], []
for r in rows:
    gid, c, a, b = r[0], r[1], int(r[2]), int(r[3])
    if gid not in want:
        continue
    E = want[gid]
    rr = runs[c]
    k_ = rr.index((a, b))
    prev = rr[k_ - 1][1] if k_ else 0
    nxt = rr[k_ + 1][0] if k_ + 1 < len(rr) else len(g[c])
    span = int((E or MAXO) + PAD)
    L, R = g[c][max(prev, a - span):a], g[c][b:min(nxt, b + span)]
    res = search(L, R, E)
    rec = dict(gap=gid, chrom=c, start=a, end=b, cause="tandem_overlap", action="join", source="", src_call="",
               fill=0, x=a, y=b, repl_L=0, repl_R=0, ins_len=0, reason="")
    if res is None:
        rec.update(action="join:not_built", reason="no_overlap_found")
        info.append([gid, E or "", "", "", "", "", len(L), len(R)])
    else:
        o, v, sv, frac = res
        info.append([gid, E or "", o, v, sv, f"{frac:.3f}", len(L), len(R)])
        if frac >= MIN_FRAC:
            rec.update(source=f"self_overlap|o={o}|ident={frac:.3f}", fill=-o, y=b + o, repl_R=o)
        else:
            rec.update(action="join:not_built", reason=f"overlap_ident_{frac:.2f}")
    if rec["action"] != "join" and E and gid in arr:
        # array over-length: the reads carry E bp less tandem array than the deposited ends; remove E bp of array,
        # first from the right end's array (next to the N run), the rest from the left end's array
        aL, aR = arr[gid]
        e = int(round(E))
        tR = min(e, aR)
        tL = e - tR
        if tL <= aL and a - tL > prev and b + tR < nxt:
            rec.update(action="join", reason="", source=f"array_trim|E={e}|L={tL}|R={tR}", fill=-e, x=a - tL,
                       y=b + tR, repl_L=tL, repl_R=tR)
            info[-1].append(f"array_trim L={tL} R={tR}")
        else:
            rec["reason"] += f";array_too_short({aL}+{aR}<{e})"
    plans.append(rec)

with open(os.path.join(out, "overlap.tsv"), "w") as fh:
    fh.write("gap\texpected_overlap\toverlap\tvotes\tvotes_smoothed\tident\tL_len\tR_len\tarray_trim\n")
    fh.write("".join("\t".join(map(str, x)) + "\n" for x in info))
# windows / junctions / candidate words exactly as repair_gaps.py plan does
by = defaultdict(list)
for p in plans:
    if p["action"] == "join":
        by[p["chrom"]].append(p)
wins, junc, cand = {}, [], set()
for c, pp in by.items():
    pp.sort(key=lambda p: p["x"])
    s, parts, prev, shift = g[c], [], 0, 0
    for p in pp:
        assert p["x"] >= prev
        parts.append(s[prev:p["x"]])
        p["nx"] = p["ny"] = p["x"] + shift
        shift -= p["y"] - p["x"]
        prev = p["y"]
    parts.append(s[prev:])
    ns = "".join(parts)
    merged = []
    for p in pp:
        v = [max(0, p["nx"] - WIN), min(len(ns), p["ny"] + WIN), [p]]
        if merged and v[0] <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], v[1]); merged[-1][2] += v[2]
        else:
            merged.append(v)
    for ws, we, pl in merged:
        wn = f"{c}:{ws}-{we}"
        wins[wn] = ns[ws:we]
        for p in pl:
            j = p["nx"]
            junc.append((p["gap"], "J", wn, j - ws))
            for t in (ns[max(0, j - NEAR):j], ns[j:j + NEAR]):
                t = t.upper()
                for i in range(0, len(t) - WK + 1, STEP):
                    k = t[i:i + WK]
                    if "N" not in k:
                        cand.add(min(k, rc(k)))
write_fa(os.path.join(out, "windows.fa"), wins)
open(os.path.join(out, "pieces.fa"), "w").close()
with open(os.path.join(out, "cand.fa"), "w") as fh:
    fh.write("".join(f">{i}\n{k}\n" for i, k in enumerate(sorted(cand))))
with open(os.path.join(out, "junctions.tsv"), "w") as fh:
    fh.write("gap\tside\twindow\tpos\n" + "".join("\t".join(map(str, j)) + "\n" for j in junc))
cols = ["gap", "chrom", "start", "end", "cause", "action", "source", "src_call", "fill", "x", "y", "repl_L", "repl_R",
        "ins_len", "reason"]
with open(os.path.join(out, "plan.tsv"), "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for p in plans:
        fh.write("\t".join(str(p[k]) for k in cols) + "\n")
print(open(os.path.join(out, "overlap.tsv")).read(), end="")
print(Counter(p["action"] for p in plans), f"{len(wins)} windows, {len(junc)} junctions, {len(cand)} words",
      file=sys.stderr)
