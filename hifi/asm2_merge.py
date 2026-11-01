#!/usr/bin/env python3
"""Apply the accepted repairs of round 1 (repair_gaps.py on dep fills / asm25) and round 2 (asm2: new assembly as fill
source, only on gaps round 1 left) to the deposited haplotype, and write the round-2 assemblies.

Both rounds are planned on the deposited coordinates (plan.tsv x/y), so the edits are combined directly; overlapping
edits are an error. Round-2 junctions were read-verified in windows that carry only round-2 edits: a round-2 edit with
a round-1 edit within WIN bp is listed (near_r1) - its crossing reads saw the deposited sequence there instead.

Usage: asm2_merge.py DEPOSITED.fa R1_DIR R2_DIR LABEL OUT_DIR
  R?_DIR: repair_gaps.py work dir of this haplotype (plan.tsv, pieces.fa, repairs.tsv); R2_DIR may be a comma list
          (e.g. r2_dep1,ov_dep1: new-source fills + overlap merges from asm2_overlap.py)
  -> OUT_DIR/LABEL.repaired2.fa (main: round-1 + round-2 insert/join), OUT_DIR/LABEL.genesearch2.fa (+ gene-search
     inserts in lower case), OUT_DIR/LABEL.repaired2.bed / .genesearch2.bed, OUT_DIR/LABEL.changes2.tsv, summary to stdout"""
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repair_gaps import read_fa, write_fa  # noqa: E402

WIN = 15000
dep, r1, r2, label, out = sys.argv[1:6]
os.makedirs(out, exist_ok=True)


def edits(d, rnd):
    pl = [l.rstrip("\n").split("\t") for l in open(os.path.join(d, "plan.tsv"))]
    ix = {k: i for i, k in enumerate(pl[0])}
    xy = {r[0]: r for r in pl[1:]}
    pieces = read_fa(os.path.join(d, "pieces.fa"))
    rp = [l.rstrip("\n").split("\t") for l in open(os.path.join(d, "repairs.tsv"))]
    jx = {k: i for i, k in enumerate(rp[0])}
    res = []
    for r in rp[1:]:
        if r[jx["accepted"]] != "yes":
            continue
        g = r[jx["gap"]]
        p = xy[g]
        res.append(dict(rnd=rnd, gap=g, chrom=p[ix["chrom"]], x=int(p[ix["x"]]), y=int(p[ix["y"]]),
                        act=r[jx["action"]], cause=r[jx["cause"]], seq=pieces.get(g, ""), src=r[jx["source"]],
                        src_call=r[jx["src_call"]], jL=r[jx["junction_reads_L"]], jR=r[jx["junction_reads_R"]]))
    return res


e1 = edits(r1, 1)
e2 = [] if r2 == "none" else sum((edits(d, 2) for d in r2.split(",")), [])   # none: round 1 only (oracle)
g1 = {e["gap"] for e in e1}
dup = [e["gap"] for e in e2 if e["gap"] in g1]
c2 = Counter(e["gap"] for e in e2)
assert all(v == 1 for v in c2.values()), "a gap repaired twice in round 2"
assert not dup, f"round 2 repeats round-1 gaps {dup[:5]}"
by = defaultdict(list)
for e in e1 + e2:
    by[e["chrom"]].append(e)
for c in by:
    by[c].sort(key=lambda e: e["x"])
    for a, b in zip(by[c], by[c][1:]):
        assert b["x"] >= a["y"], f"overlapping edits {a['gap']} {b['gap']}"
r1pos = defaultdict(list)
for e in e1:
    r1pos[e["chrom"]].append((e["x"], e["y"]))
for e in e2:
    e["near_r1"] = any(x - WIN < e["y"] and e["x"] < y + WIN for x, y in r1pos[e["chrom"]])

g = read_fa(dep)
nrun = lambda seqs: sum(len(re.findall("[Nn]{10,}", s)) for s in seqs.values())
n0, l0 = nrun(g), sum(map(len, g.values()))


def build(kind):
    res, bed = {}, []
    for c, s in g.items():
        ee = [e for e in by.get(c, []) if kind == "genesearch" or e["act"] != "genesearch"]
        if not ee:
            res[c] = s; continue
        parts, prev, shift = [], 0, 0
        for e in ee:
            ins = e["seq"].lower() if e["act"] == "genesearch" else e["seq"].upper()
            parts += [s[prev:e["x"]], ins]
            bed.append((c, e["x"] + shift, e["x"] + shift + len(ins), f"{e['gap']}|r{e['rnd']}|{e['cause']}|{e['act']}"))
            shift += len(ins) - (e["y"] - e["x"])
            prev = e["y"]
        parts.append(s[prev:])
        res[c] = "".join(parts)
    return res, bed


summ = [f"{label}: N runs >= 10 deposited {n0}, length {l0}"]
for kind, suf in (("main", "repaired2"), ("genesearch", "genesearch2")):
    seqs, bed = build(kind)
    write_fa(os.path.join(out, f"{label}.{suf}.fa"), seqs)
    with open(os.path.join(out, f"{label}.{suf}.bed"), "w") as fh:
        fh.write("".join(f"{c}\t{a}\t{b}\t{n}\n" for c, a, b, n in bed))
    summ.append(f"{label}.{suf}: N runs {nrun(seqs)}, length {sum(map(len, seqs.values()))}")
    del seqs
with open(os.path.join(out, f"{label}.changes2.tsv"), "w") as fh:
    fh.write("round\tgap\tchrom\tx\ty\taction\tcause\tins_len\tsource\tsrc_call\tjunction_reads_L\tjunction_reads_R\t"
             "near_r1\n")
    for e in sorted(e1 + e2, key=lambda e: (e["chrom"], e["x"])):
        fh.write("\t".join(map(str, [e["rnd"], e["gap"], e["chrom"], e["x"], e["y"], e["act"], e["cause"],
                                     len(e["seq"]), e["src"], e["src_call"], e["jL"], e["jR"],
                                     e.get("near_r1", "")])) + "\n")
for rnd, ee in ((1, e1), (2, e2)):
    for act in ("insert", "join", "genesearch"):
        x = [e for e in ee if e["act"] == act]
        summ.append(f"round {rnd} {act}: {len(x)} accepted, {sum(len(e['seq']) for e in x)} bp"
                    + (f", near_r1 {sum(e['near_r1'] for e in x)}" if rnd == 2 else ""))
print("\n".join(summ))
