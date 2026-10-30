#!/bin/bash
cd /home/toki/toga2/b2t/full/toga_out
head -3 rejected_items.tsv
echo "== level/reason/status counts"; awk -F'\t' 'NR>1{print $1"\t"$5"\t"$6}' rejected_items.tsv | sort | uniq -c | sort -k1,1nr | head -30
echo "== N genes -> their transcripts' loss_summary statuses"
awk -F'\t' '$1=="GENE"&&$3=="N"{print $2}' loss_summary.tsv | head -3
echo "== grep where N is written for genes/transcripts (rejection producers)"
cd /home/toki/toga2/TOGA2/src/python/modules; grep -n -E "\bN\b" results_checks.py | head; sed -n 780,820p results_checks.py
grep -rn -E "IGNORED_ITEMS|SPANNING_CLASSES" constants.py cesar_wrapper_constants.py conservation_summary.py | head
