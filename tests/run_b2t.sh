#!/bin/bash
# Real run: betulina -> trizona, TOGA2 L+M genes + 200 FI controls.
set -euo pipefail
C=/home/toki/gof/code
W=/home/toki/gof/b2t
mkdir -p $W
cd $W || exit 1
export PATH=/home/toki/toga2/env_mlc2/bin:/home/toki/miniforge3/envs/sinederella/bin:$PATH
export TMPDIR=/home/toki/tmp
T=/home/toki/toga2/b2t
TO=$T/full/toga_out
Q=/home/toki/toga2/inputs/trizona_nodot.2bit
awk -F'\t' '$1=="GENE" && ($3=="L"||$3=="M"){print $2}' $TO/loss_summary.tsv > test_genes.txt
awk -F'\t' '$1=="GENE" && $3=="FI"{print $2}' $TO/loss_summary.tsv | shuf --random-source=<(yes) -n 200 > fi_controls.txt
cat test_genes.txt fi_controls.txt > genes.txt
wc -l test_genes.txt fi_controls.txt
if [ ! -s /home/toki/gof/trizona.asm20.mmi ]; then
  twoBitToFa $Q $TMPDIR/trizona.fa
  taskset -c 80-111 minimap2 -x asm20 -t 16 -d /home/toki/gof/trizona.asm20.mmi $TMPDIR/trizona.fa 2> idx.log
  rm $TMPDIR/trizona.fa
fi
taskset -c 80-111 python3 $C/gof.py --toga $TO --ref-bed $T/sbe_ref_annotation.bed --isoforms $T/sbe_ref_isoforms.tsv \
  --ref-2bit $T/sbe.2bit --query-2bit $Q --query-index /home/toki/gof/trizona.asm20.mmi --genes genes.txt --out $W/out --threads 32
echo "== test genes (TOGA2 L/M): f1000 | f10000 | neighbours"
for c in 5 8 13; do echo "-- col $c"; awk -F'\t' -v c=$c 'NR==FNR{t[$1]=1;next} FNR>1 && ($1 in t){print $2"\t"$c}' test_genes.txt out/verdicts.tsv | sort | uniq -c | sort -k2,2 -k1,1nr; done
echo "== FI controls"
for c in 5 8 13; do echo "-- col $c"; awk -F'\t' -v c=$c 'NR==FNR{t[$1]=1;next} FNR>1 && ($1 in t){print $c}' fi_controls.txt out/verdicts.tsv | sort | uniq -c | sort -k1,1nr; done
echo done
