#!/usr/bin/env python3
"""Summaries of end_causes.tsv (betulina primary contig ends) and sites_*.tsv (trizona gap sites on betulina).
Usage: summarize.py GW_DIR TCAUSES.tsv"""
import statistics
import sys
from collections import Counter, defaultdict

d, tcause = sys.argv[1], sys.argv[2]
rows = [l.rstrip("\n").split("\t") for l in open(f"{d}/end_causes.tsv")]
h, rows = rows[0], rows[1:]
ix = {k: i for i, k in enumerate(h)}
num = lambda r, k: float(r[ix[k]]) if r[ix[k]] not in ("", "nan") else None


def show(sel, title):
    print(f"\n== {title}: {len(sel)} ends")
    by = defaultdict(list)
    for r in sel:
        by[r[ix['cause']]].append(r)
    print("cause\tn\tnear\tmid\tfar\tGC10k\thet/kb\toverhang>=3\tread_partner_end\thap1_cont>=5kb\thap2_cont>=5kb\t"
          "tri:contig\ttri:gap\ttri:seq_end\ttri:unmapped")
    for c, v in sorted(by.items(), key=lambda x: -len(x[1])):
        m = lambda k: statistics.median([x for x in (num(r, k) for r in v) if x is not None] or [float('nan')])
        sh = lambda f: f"{sum(1 for r in v if f(r)) / len(v):.2f}"
        cont = lambda k: lambda r: r[ix[k]] not in ("",) and int(r[ix[k]]) >= 5000
        tc = lambda s: lambda r: r[ix['tri_context']].startswith(s)
        print(f"{c}\t{len(v)}\t{m('near_oe'):.2f}\t{m('mid_oe'):.2f}\t{m('far_oe'):.2f}\t{m('gc10k'):.1f}\t"
              f"{m('het_per_kb'):.2f}\t{sh(lambda r: int(r[ix['n_overhang']] or 0) >= 3)}\t"
              f"{sh(lambda r: r[ix['read_partner']] != '-' and '|' in r[ix['read_partner']] and int(r[ix['n_partner']]) >= 3)}\t"
              f"{sh(cont('hap1_beyond'))}\t{sh(cont('hap2_beyond'))}\t{sh(tc('contiguous'))}\t{sh(tc('gap'))}\t"
              f"{sh(tc('seq_end'))}\t{sh(tc('unmapped'))}")
    g = Counter(r[ix['tri_context']].split(":")[-1] for r in sel if r[ix['tri_context']].startswith("gap"))
    print("trizona gap causes at these ends:", dict(g.most_common()))


big = [r for r in rows if int(r[ix['ctg_len']]) >= 1000000]
mid = [r for r in rows if 100000 <= int(r[ix['ctg_len']]) < 1000000]
small = [r for r in rows if int(r[ix['ctg_len']]) < 100000]
show(big, "ends of contigs >= 1 Mb")
show(mid, "ends of contigs 100 kb - 1 Mb")
show(small, "ends of contigs < 100 kb")
# hemizygous contigs
hemi = {r[0].split("|")[0]: int(r[ix['ctg_len']]) for r in rows if r[ix['hemi']] == "1"}
print(f"\nhemizygous contigs (median depth <= 0.65 REF, >= 200 kb): {len(hemi)}, {sum(hemi.values()) / 1e6:.1f} Mb; "
      f">= 1 Mb: {sum(1 for v in hemi.values() if v >= 1e6)} ({sum(v for v in hemi.values() if v >= 1e6) / 1e6:.1f} Mb)")
# trizona sites: status per cause class
cause = {l.split("\t")[0]: l.rstrip("\n").split("\t")[-1] for l in list(open(tcause))[1:]}
for asm in ("P", "H1", "H2"):
    st = {l.split("\t")[0]: l.split("\t")[1] for l in list(open(f"{d}/sites_{asm}.tsv"))[1:]}
    by = defaultdict(Counter)
    for s, v in st.items():
        by["random" if s.startswith("rnd") else cause.get(s, "?")][v] += 1
    print(f"\n== trizona sites on betulina {asm}: status counts")
    S = ["spanned", "split_end", "split_other", "left_only", "right_only", "none"]
    print("class\t" + "\t".join(S))
    for c, v in sorted(by.items(), key=lambda x: -sum(x[1].values())):
        print(c + "\t" + "\t".join(str(v[s]) for s in S))
