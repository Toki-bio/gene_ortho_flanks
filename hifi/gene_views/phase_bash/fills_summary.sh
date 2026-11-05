#!/bin/bash
# summary of bash_all/fills.tsv: verdicts by haplotype, source, cause, length; anomalies
set -euo pipefail
cd /home/toki/gof/phase_test/bash_all || exit 1
export LC_ALL=C
T=fills.tsv
tab() { awk -F'\t' -v c="$1" 'NR>1{print $c"\t"$NF}' $T | sort | uniq -c | awk '{print "  "$2,$3,$1}'; }
echo "== by haplotype"; tab 1
echo "== by source"; tab 9
echo "== by action"; tab 7
echo "== by cause"; tab 6
echo "== fill length by verdict (min median max)"
awk -F'\t' 'NR>1{print $NF"\t"$8}' $T | sort -k1,1 -k2,2n | awk -F'\t' '{a[$1]=a[$1]" "$2;n[$1]++} END{for(v in a){split(a[v],x," ");print "  "v,n[v],x[1],x[int((n[v]+1)/2)],x[n[v]]}}'
echo "== sums"
awk -F'\t' 'NR>1{t+=$13;o+=$14;tc+=$19;oc+=$20;ti+=$21} END{print "  T_reads",t,"O_reads",o,"T_cross",tc,"O_cross",oc,"tie_cross",ti}' $T
echo "== reads crossing per fill, target reads spanning, for tested verdicts"
awk -F'\t' 'NR>1 && ($NF=="wrong_allele"||$NF=="own_or_identical"||$NF=="undecided"){k=($22>=3?"3+":$22); print $NF"\t"k}' $T | sort | uniq -c | awk '{print "  "$2,$3,$1}'
echo "== wrong_allele: median substitutions distribution"
awk -F'\t' 'NR>1 && $NF=="wrong_allele"{m=$24; k=(m<5?"2-4":m<20?"5-19":m<100?"20-99":"100+"); print k}' $T | sort | uniq -c
echo "== wrong_allele: substitutions per kb of fill (median over fills)"
awk -F'\t' 'NR>1 && $NF=="wrong_allele"{print 1000*$24/$8}' $T | sort -n | awk '{a[NR]=$1} END{print "  n",NR,"median",a[int((NR+1)/2)],"min",a[1],"max",a[NR]}'
echo "== hap2-read check: O crossers on a fill copied from O (median x)"
awk -F'\t' 'NR>1 && ($9=="dep1"||$9=="dep2") && $27!="NA"{k=($27==0?"0":$27<2?"<2":">=2"); print k}' $T | sort | uniq -c
echo "== own_or_identical: do other-haplotype reads also show 0 (fills from the other deposited haplotype)?"
awk -F'\t' 'NR>1 && $NF=="own_or_identical"{print $9,($27=="NA"?"no_O_cross":"O_x="$27)}' $T | sort | uniq -c | sort -rn | head -8
echo "== no_target_cross: fills >=15 kb, with O crossers, with entering target reads"
awk -F'\t' 'NR>1 && $NF=="no_target_cross"{n++; if($8>=15000)l++; if($20>0)o++; if($29>0)e++; if($30!="NA"&&$30>=1)h++} END{print "  n",n,">=15kb",l+0,"O_cross>0",o+0,"T_enter>0",e+0,"T_enter x/kb>=1",h+0}' $T
echo "== examples wrong_allele, most target reads"
sort -t$'\t' -k22,22nr $T | awk -F'\t' '$NF=="wrong_allele"' | head -8 | cut -f1,2,6,7,8,9,19,20,22,23,26,27 | cut -c1-120
echo "== repaired haplotype ERROR or empty rows"
awk -F'\t' 'NR>1 && NF!=32' $T | wc -l
