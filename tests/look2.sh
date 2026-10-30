#!/bin/bash
cd /home/toki/gof/b2t
g=$(awk -F'\t' 'NR==FNR{t[$1]=1;next} FNR>1 && ($1 in t) && $8=="empty_site"{print $1}' fi_controls.txt out/verdicts.tsv)
echo "gene $g"; awk -F'\t' -v g=$g '$1==g' out/verdicts.tsv | tr '\t' '\n' | paste <(head -1 out/verdicts.tsv | tr '\t' '\n') -
grep "^$g|" out/flanks.paf | cut -f1-12
grep -P "\t$g\t|^$g\t" /home/toki/toga2/b2t/full/toga_out/orthology_classification.tsv | head -3
grep "#${g%_*}#" /home/toki/toga2/b2t/full/toga_out/query_annotation.bed | cut -f1-6,9 | head
