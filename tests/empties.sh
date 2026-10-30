#!/bin/bash
cd /home/toki/gof/b2t
awk -F'\t' 'NR==FNR{t[$1]=1;next} FNR>1 && ($1 in t) && ($5=="empty_site"||$8=="empty_site"){print $1,$2,$3,$4,$5,$6,$8,$9,$10,$13}' test_genes.txt out/verdicts.tsv
