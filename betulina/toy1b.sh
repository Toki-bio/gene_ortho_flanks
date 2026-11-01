#!/bin/bash
# Toy test of the depth-vs-GC pieces on one betulina contig and a small trizona subset.
set -euo pipefail
cd /home/toki/gof/betulina || exit 1
E=/home/toki/miniforge3/envs
export PATH=$E/sinederella/bin:$PATH
MM=$E/hifiasm/bin/minimap2
MOS=$E/lrtools/bin/mosdepth
C=/home/toki/gof/betulina/code
P=/home/toki/Sicista2026/assembly/primary_eval/hifiasm_primary.p_ctg.fa
BAM=/home/toki/Sicista2026/qc/sbet_ont.bam
VCF=/home/toki/Sicista2026/qc/calls_full/clair_auto/merge_output.vcf.gz
T1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
T2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
HIFI=/home/toki/gof/reads/ERR16751478.fastq.gz
cd toy || exit 1
# trizona: 20k HiFi reads, one chromosome pair for hap SNPs
cp $T1.fai T1.fai
(zcat $HIFI || true) | head -80000 | $MM -x map-hifi -t 16 --secondary=no $T1 - 2> mm.log | python3 $C/paf_depth.py T1.fai 10000 20 > tri.depth.bed
awk '{s+=$4} END {print "tri toy mean depth", s/NR}' tri.depth.bed
samtools faidx $T1 OZ418341.1 > h1.fa; samtools faidx $T2 OZ418081.1 > h2.fa
python3 $C/chunk_fa.py h2.fa 1000000 > h2c.fa
$MM -cx asm5 --cs -t 8 h1.fa h2c.fa 2> mm2.log | python3 $C/hap_snp_bins.py 10000 > tri.hapsnp.tsv
awk '$3>=8000 {n++; s+=$4/$3*1000} END {print "hapsnp bins", n, "mean SNP/kb", s/n}' tri.hapsnp.tsv
python3 $C/bins_gc.py h1.fa 10000 > tri.gc.bed
grep -P '^OZ418341.1\t' tri.depth.bed > tri.depth1.bed || true
python3 $C/gc_curve.py tri.depth1.bed tri.gc.bed tri.hapsnp.tsv hapsnp - 1000000 | head -3
echo TOY_DONE
