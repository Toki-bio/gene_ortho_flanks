#!/usr/bin/env python3
"""Profile genes whose 1 kb flanks were not placed (non-IG/TR).

For each gene: the betulina stretch around it that has no primary whole-genome
alignment to the query (from the TOGA2 run's aln.paf, betulina = query there),
and the neighbour-interval shrink in the query (ref gap - query gap). A deletion larger than gene + flanks shows as an
unaligned stretch ~ the shrink, with syntenic neighbours and no N run.

Usage: no_anchor_profile.py GOF_OUT TEST_GENES WGA_PAF
"""
import bisect
import collections
import re
import sys

out, test_f, paf = sys.argv[1:4]
ig = re.compile(r"^(IG[HKL]V|IG[HKL][DJ]|TR[ABGD]V|TR[ABGD]J)")
test = set(l.split()[0] for l in open(test_f) if l.strip())
rows, hdr = [], None
for l in open(out + "/verdicts.tsv"):
    f = l.rstrip("\n").split("\t")
    if hdr is None:
        hdr = f
        continue
    r = dict(zip(hdr, f))
    if r["gene"] in test and not ig.match(r["gene"]) and r["f1000_verdict"] == "no_anchor":
        rows.append(r)
contigs = set(r["ref_locus"].split(":")[0] for r in rows)

# aligned betulina intervals (primary only), merged per contig
iv = collections.defaultdict(list)
for l in open(paf):
    f = l.split("\t", 13)
    if f[0] in contigs and "tp:A:P" in l:
        iv[f[0]].append((int(f[2]), int(f[3])))
merged = {}
for c, xs in iv.items():
    xs.sort()
    m = []
    for s, e in xs:
        if m and s <= m[-1][1]:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    merged[c] = m


def unaligned_run(c, s, e):
    """Unaligned stretch around [s,e): from the end of the last aligned block
    left of the gene midpoint to the start of the first one right of it."""
    m = merged.get(c, [])
    mid = (s + e) // 2
    starts = [x[0] for x in m]
    i = bisect.bisect_right(starts, mid) - 1
    if i >= 0 and m[i][1] > mid:
        return 0, "gene_mid_aligned"
    left = m[i][1] if i >= 0 else 0
    right = m[i + 1][0] if i + 1 < len(m) else None
    return (right - left if right is not None else -1), f"{left}-{right}"


cls = collections.Counter()
print("gene\tstatus\tgene_len\tunaligned_run\trun_coords\tnb_verdict\tref_gap\tquery_gap\tshrink\tclass")
for r in rows:
    c, rest = r["ref_locus"].split(":")
    s, e = map(int, rest[:-1].split("-"))
    glen = e - s
    run, coords = unaligned_run(c, s, e)
    shrink = ""
    k = "other"
    if r["nb_query_gap"] and r["nb_ref_gap"]:
        shrink = int(r["nb_ref_gap"]) - int(r["nb_query_gap"])
    if run == 0:
        k = "gene_region_aligned_in_WGA"
    elif r["nb_verdict"].startswith("syntenic") and shrink != "" and run > glen + 2000 \
            and 0.7 <= shrink / run <= 1.3:
        k = "deletion_larger_than_flanks"
    elif r["nb_verdict"].startswith("syntenic") and shrink != "" and shrink < 0.3 * max(run, 1):
        k = "unaligned_but_interval_kept"   # diverged / rearranged, not lost
    elif not r["nb_verdict"].startswith("syntenic"):
        k = "no_syntenic_neighbours"
    cls[k] += 1
    print("\t".join(map(str, [r["gene"], r["toga_status"], glen, run, coords, r["nb_verdict"],
                              r["nb_ref_gap"], r["nb_query_gap"], shrink, k])))
print("## summary", file=sys.stderr)
for k, n in cls.most_common():
    print(f"{n}\t{k}", file=sys.stderr)
