#!/usr/bin/env python3
"""Same-animal check: at long-read heterozygous SNVs, do the Illumina reads show both alleles?
stdin: bcftools query '%CHROM\t%POS\t%REF\t%ALT\t[%AD]' from mpileup at the sites; argv[1]: sites (CHROM POS REF ALT, gz).
Counts sites with Illumina depth 10-200 (ref + long-read alt reads). Classes: both alleles (ref >= 2 and alt >= 2),
ref only, alt only. Same animal: 'both' should be the large majority (long-read false hets and low depth lower it);
an unrelated animal of the same species shares only part of its hets, so 'both' drops well below that."""
import sys, gzip
want = {}
for l in gzip.open(sys.argv[1], 'rt'):
    c, p, r, a = l.split()[:4]; want[(c, int(p))] = (r, a)
n = both = ronly = aonly = other = 0
for l in sys.stdin:
    c, p, r, alts, ad = l.rstrip('\n').split('\t')
    k = (c, int(p))
    if k not in want: continue
    wr, wa = want[k]
    if r != wr: continue
    al = [r] + alts.split(','); cnt = [int(x) if x != '.' else 0 for x in ad.split(',')]
    d = dict(zip(al, cnt)); rc = d.get(wr, 0); ac = d.get(wa, 0)
    if not 10 <= rc + ac <= 200: continue
    n += 1
    if rc >= 2 and ac >= 2: both += 1
    elif ac < 2 and rc >= 2: ronly += 1
    elif rc < 2 and ac >= 2: aonly += 1
    else: other += 1
pct = lambda x: f'{100 * x / max(1, n):.1f}%'
print(f'long-read het sites tested (Illumina depth 10-200): {n} of {len(want)}')
print(f'both alleles: {both} ({pct(both)}); ref only: {ronly} ({pct(ronly)}); alt only: {aonly} ({pct(aonly)}); other: {other}')
print('verdict:', 'consistent with the same animal' if both / max(1, n) >= 0.8 else
      ('UNCLEAR' if both / max(1, n) >= 0.6 else 'NOT the same animal (or reads/assembly mismatch) - check before polishing'))
