#!/usr/bin/env python3
"""N bases in the query interval between a gene's two neighbours (from TOGA2 one2one projections).
Usage: nb_ncheck.py PROFILE.tsv GOF_OUT TOGA_OUT QUERY_2BIT"""
import os
import subprocess
import sys

prof, out, toga, q2bit = sys.argv[1:5]
proj = {}
for l in open(os.path.join(toga, "query_annotation.bed")):
    f = l.split("\t")
    proj[f[3]] = (f[0], int(f[1]), int(f[2]))
o2o = {}
for l in list(open(os.path.join(toga, "orthology_classification.tsv")))[1:]:
    tg, tt, qg, qt, cls = l.rstrip("\n").split("\t")
    if cls == "one2one" and qt in proj and tg not in o2o:
        o2o[tg] = proj[qt]
nb = {}
hdr = None
for l in open(os.path.join(out, "verdicts.tsv")):
    f = l.rstrip("\n").split("\t")
    if hdr is None:
        hdr = f
        continue
    r = dict(zip(hdr, f))
    nb[r["gene"]] = (r["nb_left"], r["nb_right"])
seen = {}
for l in open(prof):
    f = l.rstrip("\n").split("\t")
    if f[0] == "gene":
        print(l.rstrip("\n") + "\tquery_nb_interval\tN_bases\tlongest_N_run")
        continue
    a, b = nb.get(f[0], ("", ""))
    if a not in o2o or b not in o2o or o2o[a][0] != o2o[b][0]:
        print(l.rstrip("\n") + "\t\t\t")
        continue
    c = o2o[a][0]
    s, e = min(o2o[a][2], o2o[b][2]), max(o2o[a][1], o2o[b][1])
    key = (c, s, e)
    if key not in seen:
        seq = "".join(subprocess.run(["twoBitToFa", f"{q2bit}:{c}:{s}-{e}", "stdout"], capture_output=True,
                                     text=True, check=True).stdout.splitlines()[1:]).upper()
        run = best = 0
        for ch in seq:
            run = run + 1 if ch == "N" else 0
            best = max(best, run)
        seen[key] = (seq.count("N"), best)
    print(l.rstrip("\n") + f"\t{c}:{s}-{e}\t{seen[key][0]}\t{seen[key][1]}")
