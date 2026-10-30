#!/usr/bin/env python3
"""Regions to rebuild from query (trizona) HiFi reads.

kinds:
  gap      a query N run facing an unaligned reference stretch (gene call assembly_gap,
           neighbour interval with N, profile class deletion_larger_than_flanks)
  deletion a gene called deletion by gof.py (query empty site)
  control  a query window without N at an intact one2one gene (truth = the assembly)

Writes TSV: id kind q_chr q_start q_end r_chr r_start r_end genes
(q_* = query interval to rebuild, before padding; r_* = reference interval whose
k-mers also bait reads).

Usage: make_regions.py GOF_OUT PROFILE_N.tsv QUERY_2BIT N_CONTROLS > regions.tsv
"""
import collections
import os
import random
import re
import subprocess
import sys

out, prof, q2bit, n_ctrl = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])


def fetch(tb, c, s, e):
    r = subprocess.run(["twoBitToFa", f"{tb}:{c}:{s}-{e}", "stdout"], capture_output=True, text=True, check=True)
    return "".join(r.stdout.splitlines()[1:])


rows = {}
hdr = None
for l in open(os.path.join(out, "verdicts.tsv")):
    f = l.rstrip("\n").split("\t")
    if hdr is None:
        hdr = f
        continue
    rows[f[0]] = dict(zip(hdr, f))

regions = collections.OrderedDict()
# gaps: one region per query N run inside the neighbour interval of a profiled gene
for l in open(prof):
    f = l.rstrip("\n").split("\t")
    if f[0] == "gene" or f[9] != "deletion_larger_than_flanks" or not f[12] or int(f[12]) == 0:
        continue
    g = f[0]
    rc = rows[g]["ref_locus"].split(":")[0]
    rs, re_ = map(int, f[4].split("-"))
    qc, qse = f[10].split(":")
    qs, qe = map(int, qse.split("-"))
    for m in re.finditer("N+", fetch(q2bit, qc, qs, qe).upper()):
        key = (qc, qs + m.start(), qs + m.end())
        if key in regions:
            regions[key][-1].add(g)
        else:
            regions[key] = ["gap", qc, qs + m.start(), qs + m.end(), rc, rs, re_, {g}]
# deletions called by gof.py: the query empty site (edge scale first, then 1 kb, then 10 kb)
for g, r in rows.items():
    if r.get("call") != "deletion":
        continue
    q = r.get("edge_query") or r.get("f1000_query") or r.get("f10000_query")
    qc, rest = q.split(":")
    qs, qe = map(int, rest[:-1].split("-"))
    if r.get("edge_run"):
        rc, rr = r["edge_run"].split(":")
        rs, re_ = map(int, rr.split("-"))
    else:
        rc, rr = r["ref_locus"].split(":")
        rs, re_ = map(int, rr[:-1].split("-"))
    key = (qc, qs, qe)
    if key in regions:
        regions[key][-1].add(g)
    else:
        regions[key] = ["deletion", qc, qs, qe, rc, rs, re_, {g}]
# controls: 1 kb flanks present at an FI gene, query interval without N
random.seed(3)
fi = [r for r in rows.values() if r["toga_status"] == "FI" and r["f1000_verdict"] == "present" and r["f1000_query"]]
random.shuffle(fi)
for r in fi:
    if n_ctrl == 0:
        break
    qc, rest = r["f1000_query"].split(":")
    qs, qe = map(int, rest[:-1].split("-"))
    if "N" in fetch(q2bit, qc, max(0, qs - 10000), qe + 10000).upper():
        continue
    rc, rr = r["ref_locus"].split(":")
    rs, re_ = map(int, rr[:-1].split("-"))
    regions[(qc, qs, qe)] = ["control", qc, qs, qe, rc, rs, re_, {r["gene"]}]
    n_ctrl -= 1

print("id\tkind\tq_chr\tq_start\tq_end\tr_chr\tr_start\tr_end\tgenes")
for i, v in enumerate(regions.values()):
    print(f"R{i:02d}\t" + "\t".join(map(str, v[:-1])) + "\t" + ",".join(sorted(v[-1])))
