#!/bin/bash
cd /home/toki/gof/b2t
P=/home/toki/toga2/b2t/full/aln.paf
for c in ptg000106l ptg000045l ptg000006l ptg000170l; do
  echo "== $c"
  awk -F'\t' -v c=$c '$1==c{len=$2; a+=$4-$3; t[$6]+=$4-$3} END{printf "len %d aligned-to-trizona %d (%.1f%%)\n", len, a, 100*a/len; for(k in t) print "  ", k, t[k]}' $P | sort -k2,2nr | head -5
  awk -F'\t' -v c=$c 'NR==FNR{t[$1]=1;next} FNR>1 && ($1 in t) && $3 ~ "^"c":"{print $2, $1, $5, $8}' test_genes.txt out/verdicts.tsv | awk '{s[$1]++; g=g" "$2} END{for(k in s) print "  status", k, s[k]; print "  genes:" substr(g,1,600)}'
done
