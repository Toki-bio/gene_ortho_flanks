#!/bin/bash
cd /home/toki/gof/toy
awk -F'\t' 'NR>1 && $12!="" && $12!="empty_site"{print $1, $3, $4, "run="$11, $12, "gap="$13, $14}' out_del/verdicts.tsv
echo "== truth for them"
for g in $(awk -F'\t' 'NR>1 && $12!="" && $12!="empty_site"{print $1}' out_del/verdicts.tsv); do grep -P "^$g\t" truth.tsv; done
