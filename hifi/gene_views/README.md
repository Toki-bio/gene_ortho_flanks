# Gene-by-gene package (TMEM114 first)
Steps, in the order used 2026-11-04 (all small and fast; ssearch36 on 63 kb of queries takes under a minute):
1. `make_lib.py BNK SINEBase.fa dip.fa lib.fa` : repeat library.
2. `annotate_indels.py` (contains the indel finder `runs()`) : cut every allele-specific indel >=50 bp from the alignment.
3. `ssearch36 -n -m 8 -E 1e-5 -d 10 -T 16 queries.fa lib.fa > hits.m8` (queries split to 3.5 kb pieces).
4. `parse_ssearch.py hits.m8 queries.fa` then `summarize_te.py` : tiled calls and class totals.
5. `fish_reads.sh` : locate the alleles in the HiCanu contigs, pull the reads from the indexed HiCanu BAM
   (`therioserver:/home/toki/gof/asm2/hicanu/tri.contigs.bam`, .csi), map only those reads to the deposited-genome windows.
Traps: write FASTA with LF (Windows Python writes CRLF and breaks `grep ^>name$`); coordinates must be taken from the alignment
by mapping it back, not from geneloss2/final.tsv (see docs/STATUS_2026-11-04.md section 9).
Outputs: `C:\work\Tal\sicista\genes\TMEM114\`, page https://toki-bio.github.io/Tal/sicista_genes.html.
