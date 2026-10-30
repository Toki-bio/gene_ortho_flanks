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
  Third version, ortholog replaced by 100 N: 97/100 assembly_gap, 0 empty_site.
- Real control: 200 random FI genes -> 0 empty sites at any scale
  (1 kb: 192 present, 10 kb: 194 present).

Rules learned from the tests: a 10 kb flank maps in pieces (indels), so the gene
end is placed from the piece nearest the gene; an empty site needs aligned flank
within max(100 bp, 20 % of gene) of the gene (COX14, 174 bp, was a false call
from a 1 kb extrapolation); neighbour genes test synteny only (intergenic size
change hides gene-sized deletions: 38/200 false empty sites before the change).

## Results: betulina -> trizona (therioserver:/home/toki/gof/b2t/out/verdicts.tsv)

**No confirmed deletion among the 686 TOGA2 L+M genes (2026-10-30, after the
N-run fix).** EPS8L1 and TMED4, first reported as empty sites, sit at 100-N
scaffolding gaps in trizona (284 bp and 534 bp from the site): one flank's
aligned piece ended beyond the gap and was extrapolated across it. Their
alignments in `results/b2t/aln/` show this gap; they are kept as examples of
the failure, not as deletions.

381 non-IG/TR genes, 1 kb flanks: present 185, no_anchor 108, one_flank 45,
flank_not_unique 29, assembly_gap 8, synteny_broken 4.
305 IG/TR V/D/J segments (IGH ptg000106l, IGK ptg000045l): flanks are array
copies, mostly flank_not_unique.

No-anchor genes (tests/no_anchor_profile.py + tests/nb_ncheck.py): 28 lie in
betulina stretches with no whole-genome alignment to trizona (25-285 kb) whose
size matches the shrink of the neighbour interval; 26 of them have exactly one
100-N gap in that trizona interval = assembly gaps, not deletions. Left without
an N gap: OR7E24_4190 (132 kb, olfactory receptor cluster) and TMEM114_52 (70 kb).

Rules learned: an N run of >= 10 between the anchors means assembly_gap (any
size; scaffolders write 100 N), checked between the aligned piece edges, not the
extrapolated gene ends. The toy has a third version (gene replaced by 100 N):
97/100 assembly_gap, 0 empty_site.

## Status

2026-10-30: gof.py works on the first pair; no deletion confirmed yet. Next:
anchors at the edges of a WGA-unaligned stretch, for deletions larger than the
10 kb flanks (OR7E24_4190, TMEM114_52).
