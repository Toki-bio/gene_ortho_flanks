#!/usr/bin/env python3
"""Follow the reads of the OTHER copy (copy B) along a window built from copy A (step 1c of the gap plan).

Het sites as in het_window.py (from het_window.py's cs PAF). A read alignment is copy B when it covers >= MIN_INF het
sites and carries the alternative base at >= 80 % of them; copy A at <= 20 %. For copy-B reads the script prints, per
window, B coverage per 5 kb (reads, by alignment span) and every place where a B read's alignment stops while the read
goes on (>= 1 kb more read): what the read does next -
  jump(+N)   next part of the read aligns further along the window, same strand (B lacks N bp)
  flip       next part aligns on the other strand (inversion)
  away       next part aligns to nothing in the windows (B has other sequence there)
  readend    the read itself ends (< 1 kb left)
clustered per 1 kb of window position.
Usage: phase_b.py CS_PAF WINDOWS.fa [only-window-substring]"""
import re
import sys
from collections import Counter, defaultdict

paf, win = sys.argv[1:3]
only = sys.argv[3] if len(sys.argv) > 3 else ""
MIN_DEPTH, MIN_ALT, MIN_INF, BIN = 6, 3, 3, 5000
tok = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+)")
tlen = {}
for l in open(win):
    if l.startswith(">"):
        name = l[1:].split()[0]; tlen[name] = 0
    else:
        tlen[name] += len(l.strip())
alns = []
cov = defaultdict(lambda: defaultdict(int))
subs = defaultdict(lambda: defaultdict(int))
for l in open(paf):
    f = l.rstrip("\n").split("\t")
    t = f[5]
    if only and only not in t:
        continue
    q, ql, qs, qe, st, ts, te = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], int(f[7]), int(f[8])
    cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
    alt = {}
    p = ts
    for m in tok.findall(cs):
        if m[0] == ":":
            p += int(m[1:])
        elif m[0] == "*":
            alt[p] = m[2].upper(); p += 1
        elif m[0] == "-":
            p += len(m) - 1
    alns.append((q, ql, qs, qe, st, t, ts, te, alt))
    if te - ts >= 2000:
        for b in range(ts // 1000, te // 1000 + 1):
            cov[t][b] += 1
        for x, a in alt.items():
            subs[(t, x)][a] += 1
het = defaultdict(dict)                       # t -> pos -> alt base
for (t, x), d in subs.items():
    a, n = max(d.items(), key=lambda kv: kv[1])
    c = cov[t][x // 1000]
    if c >= MIN_DEPTH and n >= MIN_ALT and 0.25 <= n / c <= 0.75:
        het[t][x] = a
hpos = {t: sorted(v) for t, v in het.items()}
import bisect
by_read = defaultdict(list)
cls = {}
for i, (q, ql, qs, qe, st, t, ts, te, alt) in enumerate(alns):
    hp = hpos.get(t, [])
    lo, hi = bisect.bisect_left(hp, ts), bisect.bisect_left(hp, te)
    inf = hp[lo:hi]
    if len(inf) >= MIN_INF:
        k = sum(1 for x in inf if alt.get(x) == het[t][x]) / len(inf)
        cls[i] = "B" if k >= 0.8 else ("A" if k <= 0.2 else "?")
    by_read[q].append(i)
for t in sorted(tlen):
    if only and only not in t:
        continue
    bcov = defaultdict(set)
    ev = defaultdict(Counter)
    for q, idx in by_read.items():
        mine = [i for i in idx if alns[i][5] == t]
        if not any(cls.get(i) == "B" for i in mine):
            continue
        ql = alns[mine[0]][1]
        segs = sorted((alns[i][2] if alns[i][4] == "+" else ql - alns[i][3],
                       alns[i][3] if alns[i][4] == "+" else ql - alns[i][2], i) for i in idx)   # all windows
        for i in mine:
            if cls.get(i) != "B":
                continue
            for b in range(alns[i][6] // BIN, alns[i][7] // BIN + 1):
                bcov[b].add(q)
        # read-forward walk: after each B segment on this window, what comes next in the read
        for j, (fs, fe, i) in enumerate(segs):
            if alns[i][5] != t or cls.get(i) != "B":
                continue
            st, ts, te = alns[i][4], alns[i][6], alns[i][7]
            for direction in ("fwd", "rev"):
                if direction == "fwd":
                    rest = ql - fe
                    pos = te if st == "+" else ts
                    nxt = [s for s in segs[j + 1:] if s[0] <= fe + 500][:1]
                else:
                    rest = fs
                    pos = ts if st == "+" else te
                    nxt = [s for s in segs[:j] if s[1] >= fs - 500][-1:]
                if pos < 2000 or pos > tlen[t] - 2000:
                    continue
                if rest < 1000:
                    what = "readend"
                elif not nxt:
                    what = "away"
                else:
                    o = alns[nxt[0][2]]
                    if o[5] != t:
                        what = f"other:{o[5].split('|')[0]}"
                    elif o[4] != st:
                        what = "flip"
                    else:
                        npos = o[6] if (direction == "fwd") == (st == "+") else o[7]
                        what = f"jump({npos - pos:+d})" if abs(npos - pos) > 1000 else "cont"
                if what != "cont":
                    ev[pos // 1000][what.split("(")[0] if what.startswith("jump") else what] += 1
    print(f"== {t}  B-reads per 5 kb: " + " ".join(str(len(bcov[b])) for b in range(tlen[t] // BIN)))
    for k in sorted(ev):
        if sum(ev[k].values()) >= 3:
            print(f"   ~{k} kb  " + " ".join(f"{w}:{n}" for w, n in ev[k].most_common()))
