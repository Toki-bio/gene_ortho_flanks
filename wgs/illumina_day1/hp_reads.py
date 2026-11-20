#!/usr/bin/env python3
"""Homopolymer tracts in annotated exons: how ambiguous is their length in the reads (the ND2 oligo-C situation)?

usage: hp_reads.py <assembly.fa> <reads.bam> <exons.bed12> <out_prefix> [minlen=5] [threads=8] [contig_filter_file]

For every maximal homopolymer run >= minlen that lies inside an exon of the BED12 (blocks), every primary read with MAPQ >= 10 that spans the run plus
2 flanking bases on each side (flank bases must match the assembly) contributes one allele: delta = (read length across the window) - (assembly length).
Per tract: n reads, counts for delta -3..+3 (others pooled), the assembly-allele fraction f0, the top non-assembly allele fraction and its delta.
Classes (set in summarise()): 'assembly-supported' f0 >= 0.8; 'ambiguous' f0 < 0.6 or a non-assembly allele >= 0.3; 'contradicted' modal allele != assembly with >= 0.5.
Frameshift-relevant: delta not a multiple of 3.
"""
import sys, collections, bisect, os
import pysam
from multiprocessing import Pool

asm, bam, bed, out = sys.argv[1:5]
minlen = int(sys.argv[5]) if len(sys.argv) > 5 else 5
nth = int(sys.argv[6]) if len(sys.argv) > 6 else 8
only = set(l.strip() for l in open(sys.argv[7])) if len(sys.argv) > 7 else None

exons = collections.defaultdict(list)
for l in open(bed):
    f = l.rstrip('\n').split('\t')
    if len(f) < 12 or (only and f[0] not in only): continue
    cs = int(f[1]); sizes = [int(x) for x in f[10].strip(',').split(',')]; starts = [int(x) for x in f[11].strip(',').split(',')]
    for sz, st in zip(sizes, starts): exons[f[0]].append((cs + st, cs + st + sz, f[3], f[8]))
FL = 2


def work(contig):
    fa = pysam.FastaFile(asm); bm = pysam.AlignmentFile(bam, 'rb')
    seq = fa.fetch(contig).upper()
    iv = sorted(set((s, e) for s, e, _, _ in exons[contig]))
    merged = []
    for s, e in iv:
        if merged and s <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], e)
        else: merged.append([s, e])
    rows = []
    for s, e in merged:
        i = s
        while i < e:
            b = seq[i]
            if b not in 'ACGT': i += 1; continue
            j = i
            while j + 1 < e and seq[j + 1] == b: j += 1
            L = j - i + 1
            if L >= minlen and i - FL >= 0 and j + FL + 1 < len(seq):
                a0, a1 = i - FL, j + FL                      # inclusive anchor positions (reference)
                counts = collections.Counter(); n = 0
                for rd in bm.fetch(contig, a0, a1 + 1):
                    if rd.is_secondary or rd.is_supplementary or rd.mapping_quality < 10 or rd.reference_start > a0 or rd.reference_end <= a1: continue
                    q0 = q1 = None
                    for qp, rp in rd.get_aligned_pairs():
                        if rp == a0: q0 = qp
                        elif rp == a1: q1 = qp; break
                    if q0 is None or q1 is None: continue
                    qs = rd.query_sequence
                    if qs[q0] != seq[a0] or qs[q1] != seq[a1]: continue
                    counts[(q1 - q0) - (a1 - a0)] += 1; n += 1
                rows.append((contig, i, b, L, n, [counts.get(d, 0) for d in range(-3, 4)], sum(v for d, v in counts.items() if d < -3 or d > 3)))
            i = j + 1
    return rows


def summarise(rows):
    res = []
    for contig, pos, b, L, n, cn, other in rows:
        if n < 10: continue
        f0 = cn[3] / n
        alt = [(cn[k] / n, k - 3) for k in range(7) if k != 3]; alt.sort(reverse=True)
        topf, topd = alt[0]
        modal = max(range(7), key=lambda k: cn[k]) - 3
        cls = 'contradicted' if (modal != 0 and cn[modal + 3] / n >= 0.5) else ('assembly-supported' if f0 >= 0.8 else ('ambiguous' if (f0 < 0.6 or topf >= 0.3) else 'weak'))
        res.append((contig, pos, b, L, n, f0, topf, topd, modal, cls))
    return res


if __name__ == '__main__':
    contigs = sorted(c for c in exons if c in pysam.FastaFile(asm).references)
    allrows = []
    with Pool(nth) as p:
        for r in p.imap_unordered(work, contigs, chunksize=2): allrows += r
    with open(out + '.tracts.tsv', 'w') as o:
        o.write('contig\tpos0\tbase\tlen\tn_reads\tdm3\tdm2\tdm1\td0\tdp1\tdp2\tdp3\tother\n')
        for contig, pos, b, L, n, cn, other in allrows: o.write(f'{contig}\t{pos}\t{b}\t{L}\t{n}\t' + '\t'.join(map(str, cn)) + f'\t{other}\n')
    res = summarise(allrows)
    with open(out + '.classes.tsv', 'w') as o:
        o.write('contig\tpos0\tbase\tlen\tn_reads\tf_assembly\ttop_alt_frac\ttop_alt_delta\tmodal_delta\tclass\n')
        for r in res: o.write('\t'.join(f'{x:.3f}' if isinstance(x, float) else str(x) for x in r) + '\n')
    lines = [f'# homopolymer tracts (>= {minlen}) in annotated exons; {len(res)} tracts with >= 10 spanning reads (of {len(allrows)} tracts)']
    lines.append('len class | base | tracts | assembly-supported | ambiguous | contradicted | weak | median f_assembly')
    for lab, lo, hi in (('5', 5, 5), ('6', 6, 6), ('7', 7, 7), ('8', 8, 8), ('9', 9, 9), ('10-11', 10, 11), ('>=12', 12, 99)):
        for bases, bl in (('AT', 'A/T'), ('CG', 'C/G')):
            sel = [r for r in res if lo <= r[3] <= hi and r[2] in bases]
            if not sel: continue
            cc = collections.Counter(r[9] for r in sel); fs = sorted(r[5] for r in sel)
            lines.append(f'{lab:6s} | {bl} | {len(sel):6d} | {cc["assembly-supported"] / len(sel):6.1%} | {cc["ambiguous"] / len(sel):6.1%} | {cc["contradicted"] / len(sel):6.1%} | {cc["weak"] / len(sel):6.1%} | {fs[len(fs) // 2]:.2f}')
    fsx = [r for r in res if r[9] in ('ambiguous', 'contradicted') and (r[7] % 3 != 0 or r[8] % 3 != 0)]
    lines.append(f'frameshift-relevant (ambiguous/contradicted with a non-multiple-of-3 allele): {len(fsx)} tracts')
    open(out + '.summary.txt', 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
