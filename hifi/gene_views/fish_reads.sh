#!/bin/bash
set -euo pipefail
cd /home/toki/gof/gene_views/TMEM114b || exit 1
rm -f reads.bam* done.flag map.log
E=/home/toki/miniforge3/envs
S=$E/sinederella/bin/samtools; M=$E/hifiasm/bin/minimap2
HC=/home/toki/gof/asm2/hicanu/tri.contigs.fasta; BAM=/home/toki/gof/asm2/hicanu/tri.contigs.bam
# 1) where do the two trizona alleles sit in the HiCanu contigs
tr -d '\015' < ../TMEM114/rows.fa | grep -A1 -E '^>TMEM114_tri_hap[12]$' | grep -v '^--' > tri_rows.fa
$M -x asm10 -c --secondary=no -t 16 $HC tri_rows.fa 2>/dev/null | awk '$10>1000{print $1,$5,$6,$8,$9,$10}' > hc_hits.txt
cat hc_hits.txt
# 2) reads already aligned there: pull them straight from the indexed BAM
: > regions.txt
awk '{print $3":"$4+1"-"$5}' hc_hits.txt > regions.txt
$S view -@ 8 -b -F 0x900 $BAM $(cat regions.txt | tr '\n' ' ') > hc_reads.bam
$S fastq -@ 8 hc_reads.bam | gzip > hc_reads.fq.gz
# 3) map just those reads to the deposited-genome slices
$M -ax map-hifi -t 16 --secondary=no ref.fa hc_reads.fq.gz 2> map.log | $S view -u -F 4 -q 20 -e 'rlen>=3000' - | $S sort -@ 8 -o reads.bam -
$S index reads.bam; $S idxstats reads.bam > idxstats.txt
echo done > done.flag
