#!/bin/bash
cd /home/toki/gof/b2t
echo "== IG/TR vs other: 1 kb | 10 kb"
for c in 5 8; do echo "-- col $c"; awk -F'\t' -v c=$c 'NR==FNR{t[$1]=1;next} FNR>1 && ($1 in t){ig=($1 ~ /^(IG[HKL]V|IG[HKL][DJ]|TR[ABGD]V|TR[ABGD]J)/)?"IG/TR":"other"; print ig"\t"$c}' test_genes.txt out/verdicts.tsv | sort | uniq -c | sort -k2,2 -k1,1nr; done
