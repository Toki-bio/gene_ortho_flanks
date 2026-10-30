#!/bin/bash
cd /home/toki/gof/toy
for v in keep del; do echo "== $v"; awk -F'\t' -v v=$v 'NR>1 && !((v=="keep"&&$5=="present"&&$8=="present")||(v=="del"&&$5=="empty_site"&&$8=="empty_site")){print $1,$3,$4,$5,$6,$7,$8,$9,$10}' out_$v/verdicts.tsv; done
echo "== truth"; cat truth.tsv | head -50 > /dev/null
for g in $(awk -F'\t' 'NR>1 && !($5=="present"&&$8=="present"){print $1}' out_keep/verdicts.tsv; awk -F'\t' 'NR>1 && !($5=="empty_site"&&$8=="empty_site"){print $1}' out_del/verdicts.tsv) ; do echo "-- $g"; grep -P "^$g\t" truth.tsv; grep "^$g|" out_keep/flanks.paf out_del/flanks.paf | cut -f1-12 | sed 's/out_//'; done | sort -u
