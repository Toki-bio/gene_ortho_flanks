#!/bin/bash
# What are the 'fill is a repeat' fills: TRF, Dfam curated consensi (blastn), self-similarity, hap1 copies, k-mer CN.
# Usage: recheck_fillrep.sh OUT_DIR [CAUSE]
set -euo pipefail
O=$1
CAUSE=${2:-fill_repeat}
mkdir -p "$O"
cd "$O" || exit 1
C=/home/toki/gof/code/hifi
G=/home/toki/gof/gf_all
DF=/home/toki/refs/dfam_current
export PATH=/home/toki/miniforge3/envs/sinederella/bin:$PATH
python3 $C/recheck_fillrep.py prep $G/all_fills.fa $G/fill_cn.tsv $G/flankcn/gap_causes.tsv $CAUSE .
if [ ! -s dfam.m8 ]; then
  if [ ! -s dfam/db.nsq ] && [ ! -s dfam/db.00.nsq ]; then
    mkdir -p dfam; makeblastdb -in $DF/curated_consensi.fa -dbtype nucl -out dfam/db > /dev/null
  fi
  blastn -task dc-megablast -query core.fa -db dfam/db -evalue 1e-5 -num_threads 16 -max_target_seqs 50 \
    -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' > dfam.m8
fi
[ -s core.trf.dat ] || trf core.fa 2 7 7 80 10 50 2000 -h -ngs > core.trf.dat || true
# self-similarity (internal duplication): minimap2 with the self diagonal skipped, PAF -> m8-like columns
minimap2 -c -x asm20 -D -P -t 8 core.fa core.fa 2> mm_self.log > self.paf
awk 'BEGIN{FS=OFS="\t"} $1==$6{print $1,$6,sprintf("%.1f",100*$10/$11),$11,$3+1,$4,$8+1,$9,0,0}' self.paf > self.m8
minimap2 -c -x asm20 -N 50 -p 0.2 -t 16 /home/toki/gof/trizona.asm20.mmi core.fa > genome.paf 2> mm.log
python3 $C/recheck_fillrep.py table . $G/all_fills.wig /home/toki/gof/ends/kmer/gc_calib.txt core.trf.dat dfam.m8 $DF/curated_meta.tsv self.m8 genome.paf > fillrep.tsv
cut -f16 fillrep.tsv | sort | uniq -c
