#!/bin/bash
# Toy for orth_sites.py (known answers) and ends_prep.py / ends_reads.py on two contigs.
set -euo pipefail
cd /home/toki/gof/betulina || exit 1
E=/home/toki/miniforge3/envs
export PATH=$E/sinederella/bin:$PATH
MM=$E/hifiasm/bin/minimap2
C=/home/toki/gof/betulina/code
P=/home/toki/Sicista2026/assembly/primary_eval/hifiasm_primary.p_ctg.fa
BAM=/home/toki/Sicista2026/qc/sbet_ont.bam
mkdir -p toy2 && cd toy2 || exit 1
samtools faidx $P ptg000002l > A.fa; samtools faidx $P ptg000005l > B.fa
cat A.fa B.fa > AB.fa; samtools faidx AB.fa
python3 $C/make_orth_toy.py A.fa B.fa > toy_flanks.fa
$MM -x asm20 -t 8 AB.fa toy_flanks.fa 2> mm.log > toy.paf
python3 $C/orth_sites.py toy.paf AB.fa.fai toy_sites.tsv
python3 $C/ends_prep.py AB.fa ab 30000 20000
wc -l ab.bins.bed ab.reads3k.bed; grep -c '>' ab.ends.fa
samtools view -q 20 -M -L ab.reads3k.bed $BAM | python3 $C/ends_reads.py $P.fai | grep -E "^end|ptg000002l|ptg000005l"
samtools bedcov -Q 20 ab.bins.bed $BAM | awk '{print $4, $5/($3-$2)}' | sort -t'|' -k2,2 -k3,3n | head -64 | tr '\n' ' '
echo; echo TOY2_DONE
