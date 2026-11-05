#!/bin/bash
# bash/awk phase test on TMED4 and the 20 control sites, compared read by read with the Python table
set -euo pipefail
E=/home/toki/miniforge3/envs
S=$E/sinederella/bin/samtools
P=/home/toki/gof/phase_test
B=$P/bash_ctl
cd "$B" || exit 1
export LC_ALL=C
cut -f2 $P/controls.reads.tsv | sort -u > names.all
wc -l names.all
/usr/local/bin/seqkit grep -j 8 -f names.all /home/toki/gof/reads/ERR16751478.fastq.gz -o ctl.fq
$S fqidx ctl.fq
cut -f1 $P/controls.reads.tsv | sort -u | xargs -P 8 -n 1 bash $B/ctl_one.sh
cat */rows_*.tsv > rows.all.tsv
awk -F'\t' -f $B/ctl_compare.awk rows.all.tsv $P/controls.reads.tsv > compare.txt
cat compare.txt
for d in TMED4 ctl_OZ418343.1_65540299; do for f in $d/sum_*.tsv; do echo "$f"; cat "$f"; done; done
rm -f ctl.fq ctl.fq.fai
echo finished
