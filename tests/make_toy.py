#!/usr/bin/env python3
"""Toy with a known answer for gof.py.

Takes N intact (FI) one2one genes, cuts their trizona region (ortholog +- PAD)
from the query genome three times: keep_<gene> (unchanged; truth = present),
del_<gene> (all query projections of the gene removed; truth = empty_site) and
gap_<gene> (gene replaced by 100 N, a scaffolding gap; truth = assembly_gap).
Writes keep.fa, del.fa, gap.fa, genes.txt into OUT.

Usage: make_toy.py TOGA_OUT REF_ISOFORMS QUERY_2BIT OUT [N]
"""
import collections
import os
import random
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from gof import twobit_sizes  # noqa: E402

PAD = 60000
toga, isoforms, q2bit, out = sys.argv[1:5]
n = int(sys.argv[5]) if len(sys.argv) > 5 else 40
os.makedirs(out, exist_ok=True)

span = {}
tx2ref = {l.split("	")[1].strip(): l.split("	")[0] for l in open(isoforms)}
for line in open(sys.argv[6] if len(sys.argv) > 6 else os.path.join(os.path.dirname(isoforms), "sbe_ref_annotation.bed")):
    f = line.split("	")
    g = tx2ref.get(f[3])
    if g:
        s0, e0 = span.get(g, (int(f[1]), int(f[2])))
        span[g] = (min(s0, int(f[1])), max(e0, int(f[2])))
status = {}
for line in open(os.path.join(toga, "loss_summary.tsv")):
    f = line.rstrip("\n").split("\t")
    if f[0] == "GENE":
        status[f[1]] = f[2]
proj = {}
for line in open(os.path.join(toga, "query_annotation.bed")):
    f = line.rstrip("\n").split("\t")
    proj[f[3]] = (f[0], int(f[1]), int(f[2]))
tx2g, o2o = {}, {}
lines = open(os.path.join(toga, "orthology_classification.tsv")).read().splitlines()[1:]
for line in lines:
    tg, tt, qg, qt, cls = line.split("\t")
    tx2g[tt] = tg
    if cls == "one2one" and qt in proj:
        o2o.setdefault(tg, proj[qt])
gproj = collections.defaultdict(list)
for name, loc in proj.items():
    tx = "#".join(name.split("#")[:2])
    if tx in tx2g:
        gproj[tx2g[tx]].append(loc)

sizes = twobit_sizes(q2bit)

cand = []
for g, (c, s, e) in o2o.items():
    if status.get(g) != "FI":
        continue
    locs = [p for p in gproj[g] if p[0] == c and p[1] < e + PAD and p[2] > s - PAD]
    gs, ge = min(p[1] for p in locs), max(p[2] for p in locs)
    if not (5000 <= ge - gs <= 100000) or gs - PAD < 0 or ge + PAD > sizes[c]:
        continue
    # truth needs the cut to remove the whole reference gene: query span within 0.8-1.25 of it
    rlen = span[g][1] - span[g][0]
    if not 0.8 <= (ge - gs) / rlen <= 1.25:
        continue
    cand.append((g, c, gs, ge))
random.seed(1)
random.shuffle(cand)
# toy regions must not overlap, or a flank sits on two toy contigs (MAPQ 0)
pick, taken = [], []
for g, c, gs, ge in cand:
    if any(c == c2 and gs - PAD < e2 and ge + PAD > s2 for c2, s2, e2 in taken):
        continue
    pick.append((g, c, gs, ge))
    taken.append((c, gs - PAD, ge + PAD))
    if len(pick) == n:
        break
pick.sort()


def fetch(c, s, e):
    r = subprocess.run(["twoBitToFa", f"{q2bit}:{c}:{s}-{e}", "stdout"], capture_output=True, text=True, check=True)
    return "".join(r.stdout.splitlines()[1:])


with open(os.path.join(out, "keep.fa"), "w") as k, open(os.path.join(out, "del.fa"), "w") as d, \
        open(os.path.join(out, "gap.fa"), "w") as gp, open(os.path.join(out, "genes.txt"), "w") as gl, open(os.path.join(out, "truth.tsv"), "w") as t:
    t.write("gene\tquery_locus\tdeleted_len\n")
    for g, c, gs, ge in pick:
        left, mid, right = fetch(c, gs - PAD, gs), fetch(c, gs, ge), fetch(c, ge, ge + PAD)
        k.write(f">keep_{g}\n{left}{mid}{right}\n")
        d.write(f">del_{g}\n{left}{right}\n")
        gp.write(f">gap_{g}\n{left}{'N' * 100}{right}\n")
        gl.write(g + "\n")
        t.write(f"{g}\t{c}:{gs}-{ge}\t{ge - gs}\n")
print(f"{len(pick)} toy genes of {len(cand)} candidates -> {out}")
