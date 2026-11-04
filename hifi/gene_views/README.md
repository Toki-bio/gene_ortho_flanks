# Gene-by-gene package (TMEM114 first)
1. `make_lib.py BNK SINEBase.fa dip.fa lib.fa`: repeat library.
2. Indels: pairwise `minimap2 -x asm20 -c --cs=short -r 2000,120000 -z 1000000` between the sequences, every + or - of 50 bp or more
   (code in audit_2026-11-04/a4.sh, query list in a5.sh). Do not take indels from MAFFT columns: it forces unrelated ends together,
   and a sequence that simply ends looks like an indel.
3. `ssearch36 -n -m 8 -E 1e-5 -T 16 q.fa lib.fa > q.m8` (seconds), then `parse_ssearch.py` (needs lib.fa, q.m8 and the .bnk for species).
4. `fish_reads.sh`: locate alleles in HiCanu contigs, pull reads from the indexed BAM, map only those to the deposited windows.

Traps:
- geneloss2/final.tsv hapN_locus = exons +-1000 bp in the REPAIRED assemblies (repair/depN.repaired.fa), 0-based; not deposited coordinates.
- PAF starts are 0-based: region_start + paf_start is already the 1-based genomic position.
- het_genes/*.bed are per-sequence ungapped positions, not alignment columns. Check: exon 1 must start ATG.
- trizona.asm20.mmi names sequences without the .1 version.
- Windows Python writes CRLF unless newline="\n".

Outputs: C:\work\Tal\sicista\genes\TMEM114\, page https://toki-bio.github.io/Tal/sicista_genes.html.
