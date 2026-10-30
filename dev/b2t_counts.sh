#!/bin/bash
set -euo pipefail
cd /home/toki/toga2/b2t/full/toga_out
echo "== gene-level status"; awk -F'\t' '$1=="GENE"{c[$3]++} END{for(k in c) print k,c[k]}' loss_summary.tsv | sort
echo "== genes by orthology class (t_gene unique)"; awk -F'\t' 'NR>1{print $1"\t"$5}' orthology_classification.tsv | sort -u | cut -f2 | sort | uniq -c
echo "== query_genes.bed"; head -2 query_genes.bed
