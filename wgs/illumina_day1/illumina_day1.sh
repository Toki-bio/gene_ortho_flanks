#!/bin/bash
# Day-one checks and polishing for the S. betulina Sb1 Illumina reads (prepared 2026-11-21, before the reads arrived).
# usage: illumina_day1.sh R1.fq.gz R2.fq.gz OUT_DIR [threads=48] [asm=primary] [het_sites.tsv.gz] [tract_classes.tsv] [exons.bed12]
# Steps (each skipped if its output exists, so the script can be rerun after a failure):
#  1 fastp                      -> trim_1/2.fq.gz, fastp.json
#  2 bwa-mem2 + fixmate/markdup -> ill.bam(.bai), flagstat.txt, depth.txt (median depth over the assembly)
#  3 same-animal check          -> same_animal.txt (long-read het sites: are both alleles in the Illumina reads?)
#  4 Illumina k-mer QV          -> qv_before.txt (meryl k=21; Merqury formula)
#  5 homopolymer tracts         -> ill_tracts.* (hp_reads.py on the Illumina BAM) + tract_compare.txt (vs long-read classes)
#  6 pypolca --careful          -> polished.fa, pypolca report; then QV again -> qv_after.txt
set -euo pipefail
R1=$(readlink -f $1); R2=$(readlink -f $2); O=$3; T=${4:-48}
A=$(readlink -f ${5:-/home/toki/Sicista2026/illumina_prep/index/primary.fa})
SITES=$(readlink -f ${6:-/home/toki/Sicista2026/illumina_prep/het_sites_500k.tsv.gz})
CLS=$(readlink -f ${7:-/home/toki/Sicista2026/demog/hp/bet.classes.tsv})
BED=$(readlink -f ${8:-/home/toki/Sicista2026/hiller_share/TOGA2/query_annotation.bed})
HERE=$(dirname $(readlink -f $0))
mkdir -p $O; cd $O || exit 1
export PATH=/home/toki/miniforge3/envs/polish/bin:$PATH TMPDIR=/home/toki/tmp
PY=/home/toki/Sicista2026/demog/envs/smcpp/bin/python          # has pysam
step() { echo "#### $* -- $(date '+%F %T')"; }
[ -s $A.bwt.2bit.64 ] || { step "bwa-mem2 index (missing)"; bwa-mem2 index $A > index.log 2>&1; }

step "1 fastp"
[ -s fastp.json ] || fastp -i $R1 -I $R2 -o trim_1.fq.gz -O trim_2.fq.gz -w 16 -j fastp.json -h fastp.html 2> fastp.log

step "2 map"
if [ ! -s ill.bam.bai ]; then
  bwa-mem2 mem -t $T -R '@RG\tID:sb1ill\tSM:sb1' $A trim_1.fq.gz trim_2.fq.gz 2> bwa.log \
    | samtools fixmate -m -@4 - - | samtools sort -@8 -m 2G -T $TMPDIR/ill_sort - \
    | samtools markdup -@8 - ill.bam 2> markdup.log
  samtools index -@8 ill.bam
fi
[ -s flagstat.txt ] || samtools flagstat -@8 ill.bam > flagstat.txt
[ -s depth.txt ] || samtools depth -a -@4 ill.bam | awk 'NR%100==1{print $3}' | sort -n | awk '{a[NR]=$1} END{print "median depth (1/100 positions):", a[int(NR/2)+1]}' > depth.txt

step "3 same-animal"
[ -s same_animal.txt ] || { bcftools mpileup -f $A -T $SITES -a FORMAT/AD -q 20 -Q 20 -d 500 --threads 8 ill.bam 2>/dev/null \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' | $PY $HERE/same_animal.py $SITES > same_animal.txt; }

step "4 QV before"
[ -s qv_before.txt ] || bash $HERE/qv_meryl.sh $A trim_1.fq.gz trim_2.fq.gz before $T > qv_before.txt

step "5 tracts"
[ -s ill_tracts.classes.tsv ] || $PY $HERE/hp_reads.py $A ill.bam $BED ill_tracts 5 $T > hp_reads.log 2>&1
[ -s tract_compare.txt ] || $PY $HERE/tract_compare.py $CLS ill_tracts.classes.tsv > tract_compare.txt

step "6 pypolca"
[ -s polished.fa ] || { pypolca run -a $A -1 trim_1.fq.gz -2 trim_2.fq.gz -t $T -o pypolca --careful -f > pypolca.log 2>&1; cp pypolca/pypolca_corrected.fasta polished.fa; }
[ -s qv_after.txt ] || bash $HERE/qv_meryl.sh polished.fa trim_1.fq.gz trim_2.fq.gz after $T > qv_after.txt
step "done"
for f in flagstat.txt depth.txt same_animal.txt qv_before.txt tract_compare.txt qv_after.txt; do echo "== $f"; cat $f; done
