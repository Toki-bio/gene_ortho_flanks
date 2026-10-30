#!/bin/bash
cd /home/toki/gof/b2t
awk -F'\t' 'NR==FNR{t[$1]=1;next} FNR>1 && ($1 in t){ig=($1 ~ /^(IG[HKL]V|IG[HKL][DJ]|TR[ABGD]V|TR[ABGD]J)/)?"IG/TR_segment":"other"; print ig"\t"$5}' test_genes.txt out/verdicts.tsv | sort | uniq -c | sort -k2,2 -k1,1nr
echo "== other genes, 1 kb verdict x neighbour verdict"
awk -F'\t' 'NR==FNR{t[$1]=1;next} FNR>1 && ($1 in t) && $1 !~ /^(IG[HKL]V|IG[HKL][DJ]|TR[ABGD]V|TR[ABGD]J)/{print $5"\t"$13}' test_genes.txt out/verdicts.tsv | sort | uniq -c | sort -k1,1nr | head -20
