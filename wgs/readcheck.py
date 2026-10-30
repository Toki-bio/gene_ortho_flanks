#!/usr/bin/env python3
"""Read support for gapfill.py fills: reads mapped (PAF, map-hifi) to each haplotype assembly
separately; for a fill contig:lo-hi a read supports the left / right junction when one
alignment with MAPQ >= MIN_MAPQ and no large indel (read span ~ contig span) covers lo-M..lo+M / hi-M..hi+M, and spans the fill when it
covers lo-M..hi+M. Adds to gaps.tsv: <hap>_reads_left, _reads_right, _reads_span, and
<hap>_supported (both junctions >= MIN_READS reads).

Usage: readcheck.py OUT_DIR(of gapfill.py) LABEL=reads_vs_hap.paf [LABEL=...]"""
import os
import sys

M = 500
MIN_MAPQ = 5
MIN_READS = 3
out = sys.argv[1]
pafs = dict(a.split("=", 1) for a in sys.argv[2:])

rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(out, "gaps.tsv"))]
hdr, rows = rows[0], rows[1:]
for label, paf in pafs.items():
    ci = hdr.index(f"{label}_contig")
    si = hdr.index(f"{label}_status")
    want = {}                                    # contig -> [(lo, hi, row index)]
    for i, r in enumerate(rows):
        if r[si] == "spanned":
            t, se = r[ci].rsplit(":", 1)
            lo, hi = map(int, se.split("-"))
            want.setdefault(t, []).append((lo, hi, i))
    cnt = {i: [0, 0, 0] for v in want.values() for _, _, i in v}
    for l in open(paf):
        f = l.split("\t", 12)
        t = f[5]
        if t not in want or int(f[11]) < MIN_MAPQ:
            continue
        ts, te = int(f[7]), int(f[8])
        # a read of the other haplotype can cross a junction inside one alignment with a large
        # indel (het toy: a 5 kb deletion doubled the count); keep alignments without one
        if abs((te - ts) - (int(f[3]) - int(f[2]))) > max(200, 0.02 * (te - ts)):
            continue
        for lo, hi, i in want[t]:
            left = ts <= lo - M and te >= lo + M
            right = ts <= hi - M and te >= hi + M
            c = cnt[i]
            c[0] += left
            c[1] += right
            c[2] += ts <= lo - M and te >= hi + M
    hdr += [f"{label}_reads_left", f"{label}_reads_right", f"{label}_reads_span", f"{label}_supported"]
    for i, r in enumerate(rows):
        if i in cnt:
            a, b, s = cnt[i]
            r += [a, b, s, "yes" if min(a, b) >= MIN_READS else "no"]
        else:
            r += ["", "", "", ""]
with open(os.path.join(out, "gaps_reads.tsv"), "w") as h:
    for r in [hdr] + rows:
        h.write("\t".join(map(str, r)) + "\n")
