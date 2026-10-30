#!/bin/bash
# L/M genes with an empty site at any flank scale (1 kb, 10 kb, edge of the WGA-unaligned stretch)
cd /home/toki/gof/b2t || exit 1
awk -F'\t' 'NR==FNR{t[$1]=1;next} FNR>1 && ($1 in t) && ($5=="empty_site"||$8=="empty_site"||$12=="empty_site"){print $1,$2,$3,$4,"f1k="$5,"f10k="$8,"edge="$12,"run="$11,"gap="$13,$14,"nb="$17}' test_genes.txt out/verdicts.tsv
