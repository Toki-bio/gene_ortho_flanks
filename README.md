# gene_ortho_flanks

Flank-based orthology and gene-loss evidence to supplement TOGA2, built on the
same logic as [SINE_orth_loc](https://github.com/Toki-bio/SINE_orth_loc): a locus
is placed in the other genome by its flanks, not by itself. An empty site (flanks
adjacent, gene absent) is direct evidence of a deletion.

Intended for close genome pairs (non-coding flanks must align). Not for
human-vs-rodent TOGA2 runs.

## Inputs (first pair: *Sicista betulina* -> *S. trizona*)

- TOGA2 run: `therioserver:/home/toki/toga2/b2t/full/toga_out/` (v2.0.10a;
  reference = betulina hifiasm primary, annotation = intact betulina projections
  from Hiller's TOGA2 vs hg38, `../../sbe_ref_annotation.bed`; query = trizona).
- Gene-level status (loss_summary.tsv): FI 13,002, I 1,422, PI 975, UL 369,
  L 112, M 574, N 1,247, PG 42.
- Orthology (genes): one2one 15,158, one2zero 1,975, many2one 309,
  many2many 193, one2many 108.

## Method

Two anchor scales per reference gene, left and right:

1. **Neighbour genes**: nearest flanking reference genes with a one2one FI/I
   ortholog in the query (TOGA2's own calls); gives the expected query interval.
2. **Sequence flanks**: the nearest flank (next to the gene) and the flank 10 kb
   out, mapped to the query with minimap2 (settings tested in SINE_orth_loc).

Per gene the verdict is the following:

| anchors in query | verdict |
|---|---|
| same contig, order/strand right, gap = reference span minus gene | deletion confirmed (empty site) |
| same contig, gap = full reference span | locus present; L/UL call is decay, check mutations with reads |
| an anchor at a contig end or N-run | unresolved: assembly |
| different contigs / wrong order | synteny broken, no call |

For one2many/many2many: the query copy inside the anchor interval is the
ortholog. Only one partner is allowed per query site, because SINE_orth_loc's
rescue step paired both copies of a duplication with one site.

Each confirmed deletion gets a reference-locus+flanks vs query-site alignment
for viewing.

## TOGA2 statuses used (from TOGA2 v2.0.10 source)

- **N** is the lowest rank (-1). The reference transcript was rejected at setup,
  before any projection: CDS length not a multiple of 3 (`OUT_OF_FRAME`,
  input_producer.py ~l.839; 3,065 transcripts here) or `ILLEGAL_NAME` (202).
  It says nothing about the query, so N genes are **not** deletion candidates.
- **M** at transcript level: `NO_PROJ` (no classifiable projection, 527) or
  `ALL_ORTHS_DISCARDED` (1,092). These are the main deletion candidates.
- **L** projections with `NO_EXONS_ALIGNED` (96): a chain exists but no exon
  aligned. These are candidates for an empty site inside a syntenic chain.

Test set: gene-level L (112) + M (574).

## Tests

- Toy with a known answer (`tests/make_toy.py`, `tests/run_toy.sh`): 100 intact
  one2one genes, trizona region +-60 kb kept or with the ortholog cut out.
  97/100 right at both flank scales, the other 3 refused, no wrong calls.
- Real control: 200 random FI genes -> 0 empty sites at any scale
  (1 kb: 192 present, 10 kb: 194 present).

Rules learned from the tests: a 10 kb flank maps in pieces (indels), so the gene
end is placed from the piece nearest the gene; an empty site needs aligned flank
within max(100 bp, 20 % of gene) of the gene (COX14, 174 bp, was a false call
from a 1 kb extrapolation); neighbour genes test synteny only (intergenic size
change hides gene-sized deletions: 38/200 false empty sites before the change).

## Results: betulina -> trizona (therioserver:/home/toki/gof/b2t/out/verdicts.tsv)

686 TOGA2 L+M genes, 1 kb flanks: present 211, no_anchor 206, one_flank 141,
size_changed 69, synteny_broken 54, empty_site 2. Neighbours syntenic 504.
Empty sites, both flank scales agree and neighbours syntenic:
EPS8L1 (L, 17,952 bp gene -> 272 bp in trizona), TMED4 (M, 2,947 bp -> 299 bp).

Alignments (ViewAlign links): `results/b2t/aln/`. EPS8L1: one 17,679-column gap
in trizona, 220 of 2,173 CDS bases still opposite trizona bases (3' end), flank
identity 95.8 / 97.8 %. TMED4: one 2,652-column gap, 150 of 684 CDS bases
opposite trizona bases, flanks 97.9 %. Both are deletions of most of the gene.

Why flanks are not placed (tests/why_unplaced.py, tests/split_ig.sh): 305 of the
686 L/M genes are IG/TR V/D/J segments (IGH array on ptg000106l, IGK on
ptg000045l). Their flanks are array copies (MAPQ < 20 or a paralogous copy), so
size_changed / synteny_broken there are not trustworthy. The other 381 genes:
present 194, no_anchor 108 (61 of them with neighbours syntenic but the
neighbour interval changed size), one_flank 62, empty_site 2.

## Status

2026-10-30: gof.py works on the first pair; next: why L/M flanks fail to place
(no_anchor / one_flank), and alignments for the empty sites.
