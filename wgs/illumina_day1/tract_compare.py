#!/usr/bin/env python3
"""Exon homopolymer tracts: long-read class (hp/bet.classes.tsv from the CycloneSEQ reads) against the Illumina class
(hp_reads.py on the Illumina BAM). usage: tract_compare.py long_classes.tsv illumina_classes.tsv
Illumina verdict per tract (>= 10 spanning Illumina reads): 'assembly right' f_assembly >= 0.8; 'assembly wrong' modal allele
!= assembly with >= 0.8 of reads (pypolca should fix these); 'Illumina ambiguous' otherwise; 'no Illumina call' < 10 reads.
Main number: of the long-read ambiguous/weak/contradicted tracts of 8 bp or less, how many Illumina settles."""
import sys, collections
def load(p):
    d = {}
    for l in open(p):
        f = l.rstrip('\n').split('\t')
        if f[0] == 'contig': h = f; continue
        d[(f[0], f[1])] = dict(zip(h, f))
    return d
lr, il = load(sys.argv[1]), load(sys.argv[2])
def verdict(t):
    if t is None or int(t['n_reads']) < 10: return 'no Illumina call'
    fa = float(t['f_assembly']); md = int(t['modal_delta'])
    if fa >= 0.8: return 'assembly right'
    if md != 0 and 1 - fa >= 0.8 and float(t['top_alt_frac']) >= 0.8: return 'assembly wrong'
    return 'Illumina ambiguous'
tab = collections.Counter(); fs = collections.Counter()
for k, t in lr.items():
    L = int(t['len']); lc = t['class']; band = '<=8 bp' if L <= 8 else '>=9 bp'
    v = verdict(il.get(k)); tab[(lc, band, v)] += 1
    if v == 'assembly wrong':
        d = int(il[k]['modal_delta']); fs['frameshift' if d % 3 else 'in-frame'] += 1
V = ['assembly right', 'assembly wrong', 'Illumina ambiguous', 'no Illumina call']
print('long-read class | length | ' + ' | '.join(V))
for lc in ('assembly-supported', 'weak', 'ambiguous', 'contradicted'):
    for band in ('<=8 bp', '>=9 bp'):
        row = [tab[(lc, band, v)] for v in V]
        if sum(row): print(f'{lc} | {band} | ' + ' | '.join(map(str, row)))
unsettled = sum(tab[(lc, '<=8 bp', v)] for lc in ('weak', 'ambiguous', 'contradicted') for v in V)
settled = sum(tab[(lc, '<=8 bp', v)] for lc in ('weak', 'ambiguous', 'contradicted') for v in V[:2])
print(f'long-read unsettled tracts <= 8 bp: {unsettled}; settled by Illumina: {settled} ({100 * settled / max(1, unsettled):.1f}%)')
print(f'tracts where Illumina says the assembly is wrong: {sum(fs.values())} (frameshift {fs["frameshift"]}, in-frame {fs["in-frame"]})')
