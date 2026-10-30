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

## Status

2026-10-30: repository started; inputs located.
