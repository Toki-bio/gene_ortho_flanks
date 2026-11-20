# Day-one pipeline for the S. betulina Sb1 Illumina reads

Prepared 2026-11-21, before the reads arrived. Server copy: `therioserver:/home/toki/Sicista2026/illumina_prep/day1/`.
Ready inputs on therioserver:
- bwa-mem2 index of the betulina primary assembly: `/home/toki/Sicista2026/illumina_prep/index/primary.fa.*`
  (symlink to `assembly/primary_eval/hifiasm_primary.p_ctg.fa`; built 2026-11-21, log `index.log`, marker `INDEX_DONE`).
- 500,467 random long-read het SNVs (PSMC set: PASS, 0/1, AF 0.3-0.7): `/home/toki/Sicista2026/illumina_prep/het_sites_500k.tsv.gz`.
- Long-read tract classes: `demog/hp/bet.classes.tsv` (built with Hiller's `hiller_share/TOGA2/query_annotation.bed`; use the same BED).

Run, when the reads arrive:

    bash illumina_day1.sh R1.fq.gz R2.fq.gz OUT_DIR 48

| step | output | what it answers |
|---|---|---|
| 1 fastp | `fastp.json` | read quality |
| 2 bwa-mem2, fixmate, sort, markdup | `ill.bam`, `flagstat.txt`, `depth.txt` | mapping rate, depth |
| 3 same-animal | `same_animal.txt` | are both alleles of the long-read het sites in the Illumina reads? |
| 4 QV before | `qv_before.txt` | assembly QV from Illumina k-mers (meryl k=21, Merqury formula) |
| 5 tracts | `ill_tracts.*`, `tract_compare.txt` | which long-read-ambiguous exon homopolymers the Illumina reads settle, and where the assembly is wrong |
| 6 pypolca --careful, QV after | `polished.fa`, `qv_after.txt` | polished primary |

Every step is skipped if its output exists, so a failed run can be restarted.

## Toy test (`toy_day1.sh`, contig ptg000170l, 10 Mb)
Truth: hap A = contig, hap B = contig + its 23,846 long-read het SNVs; 30x 2x150 wgsim reads; 20 planted 1-bp deletions in exon
homopolymers of the toy assembly. Result (2026-11-21):
- same-animal: 96.4% of sites show both alleles (verdict threshold 80%);
- QV 58.18 before, 64.31 after polishing;
- pypolca restored 19 of 20 planted errors;
- tract check: all 6 planted tracts still >= 5 bp called 'assembly wrong' (frameshift); 14 shortened to 4 bp fall below the 5-bp scan floor
  (they are still fixed by pypolca, but not listed by the tract check); 98% of long-read tracts get an Illumina call.
Bugs found by the toy and fixed: meryl statistics label (`present`), exon BED must be Hiller's (the other one covers ~45% of the tracts).
