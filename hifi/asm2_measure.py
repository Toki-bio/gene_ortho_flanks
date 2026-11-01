#!/usr/bin/env python3
"""Round-2 gap census: for every gap of one deposited haplotype that repair round 1 left unrepaired (repairs.tsv
accepted != yes), does a new assembly span both flanks (wgs/gapfill.py statuses), and what is the fill's GC-calibrated
k-mer copy-number call (fill_cn.py)? Also writes the cause file for repair round 2 (asm2_repair2.sh):
  round-1 causes kept as they are (hom_one_hap, fill_nocall, join, one_copy_fill, fill_unclear, fill_short_or_unscored,
  fill_repeat) - the same rule then decides on the new fill;
  'nobody spans' causes (tandem, high_copy, other, ga_tract) -> fill_nocall (the new fill's copy-number call decides:
  two_copies -> insert, all repeat -> leave, else gene-search only);
  chrX -> chrX_fill (male X has one copy: a non-repeat fill is the only copy and goes into the main assembly).

Usage: asm2_measure.py GAPS.tsv FILL_CN.tsv FILLS.fa CAUSES.tsv REPAIRS.tsv ASSEMBLY OUT_PREFIX
  -> OUT_PREFIX.remaining.tsv (one row per unrepaired gap), OUT_PREFIX.by_cause.tsv, OUT_PREFIX.causes2.tsv,
     OUT_PREFIX.gaps2.txt (comma list of unrepaired gaps with >= 1 new spanning source)"""
import sys
from collections import Counter, defaultdict

gtsv, cntsv, ffa, causes, repairs, asm, out = sys.argv[1:8]
NOBODY = {"tandem", "high_copy", "other", "ga_tract"}

rows = [l.rstrip("\n").split("\t") for l in open(gtsv)]
h, rows = rows[0], rows[1:]
ix = {k: i for i, k in enumerate(h)}
labels = [c[:-7] for c in h if c.endswith("_status")]
heads = [l[1:].strip() for l in open(ffa) if l.startswith(">")]
mine = set(heads)
call = {}
cnh = None
for l in open(cntsv):
    f = l.rstrip("\n").split("\t")
    if cnh is None:
        cnh = f; continue
    call.setdefault((f[0], f[1]), []).append(f[6])
# FILL_CN.tsv must be computed on this assembly's FILLS.fa only (gap ids repeat between dep1 and dep2)
assert sum(len(v) for v in call.values()) == len(heads), "FILL_CN.tsv does not match FILLS.fa"
cause_h = None
cause = {}
for l in open(causes):
    f = l.rstrip("\n").split("\t")
    if cause_h is None:
        cause_h = f; continue
    cause[f[0]] = f[-1]
acc = {}
for l in list(open(repairs))[1:]:
    f = l.rstrip("\n").split("\t")
    if f[0] == asm:
        acc[f[1]] = (f[6], f[13], f[14])          # action, accepted, reason


def cn_of(gid, lab):
    v = call.get((gid, lab), [])
    if not v:
        return ""
    return v[0] if len(v) == 1 else "|".join(sorted(set(v)))


rem = []
for r in rows:
    gid = r[0]
    a = acc.get(gid, ("", "no", "not in round 1"))
    if a[1] == "yes":
        continue
    st = {l: r[ix[f"{l}_status"]] for l in labels}
    fl = {l: r[ix[f"{l}_fill"]] for l in labels}
    rem.append(dict(gap=gid, chrom=r[1], cause=cause.get(gid, "other"), r1_action=a[0], r1_reason=a[2], st=st, fl=fl,
                    cn={l: (cn_of(gid, l) if st[l] == "spanned" else "") for l in labels}))

with open(out + ".remaining.tsv", "w") as fh:
    fh.write("\t".join(["gap", "chrom", "cause", "r1_action", "r1_reason"] +
                       sum([[f"{l}_status", f"{l}_fill", f"{l}_cn"] for l in labels], [])) + "\n")
    for x in rem:
        fh.write("\t".join([x["gap"], x["chrom"], x["cause"], x["r1_action"], x["r1_reason"]] +
                           sum([[x["st"][l], x["fl"][l], x["cn"][l]] for l in labels], [])) + "\n")

by = defaultdict(list)
for x in rem:
    by[x["cause"]].append(x)
with open(out + ".by_cause.tsv", "w") as fh:
    fh.write("cause\tremaining\t" + "\t".join(f"{l}_spanned" for l in labels) + "\tany_spanned\t" +
             "\t".join(f"{l}_calls" for l in labels) + "\n")
    tot = Counter()
    for cs in sorted(by):
        xx = by[cs]
        sp = {l: sum(x["st"][l] == "spanned" for x in xx) for l in labels}
        anys = sum(any(x["st"][l] == "spanned" for l in labels) for x in xx)
        calls = {l: Counter(x["cn"][l] for x in xx if x["st"][l] == "spanned") for l in labels}
        fh.write(f"{cs}\t{len(xx)}\t" + "\t".join(str(sp[l]) for l in labels) + f"\t{anys}\t" +
                 "\t".join(",".join(f"{k}:{v}" for k, v in sorted(calls[l].items())) for l in labels) + "\n")
        tot["n"] += len(xx); tot["any"] += anys
        for l in labels:
            tot[l] += sp[l]
    fh.write(f"TOTAL\t{tot['n']}\t" + "\t".join(str(tot[l]) for l in labels) + f"\t{tot['any']}\n")

with open(out + ".causes2.tsv", "w") as fh:
    fh.write("gap\tcause\n")
    for r in rows:
        c = cause.get(r[0], "other")
        if c in NOBODY:
            c = "fill_nocall"
        elif c == "chrX":
            c = "chrX_fill"
        fh.write(f"{r[0]}\t{c}\n")
g2 = [x["gap"] for x in rem if any(x["st"][l] == "spanned" for l in labels)]
open(out + ".gaps2.txt", "w").write(",".join(g2) + "\n")
print(open(out + ".by_cause.tsv").read(), end="")
print(f"{len(g2)} unrepaired gaps with a new spanning source", file=sys.stderr)
