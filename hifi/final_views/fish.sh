#!/bin/bash
# one fishing run over the whole HiFi file for all windows of both haplotypes, reads then allotted to items
set -euo pipefail
. "$(dirname "$0")/cfg.sh"
cd "$W" || exit 1
for d in w/*; do
	. $d/meta
	# T near the gap, O whole, F fill
	awk -v a=$g0 -v b=$g1 '!/^>/{ l=a<15000?a:15000; print substr($0,a-l+1,l); print substr($0,b+1,15000) }' $d/T.fa
	if [ -s $d/O.fa ]; then grep -v '>' $d/O.fa; fi
	if [ "$closed" = 1 ]; then awk -v a=$f0 -v b=$f1 '!/^>/&&b>a{print substr($0,a+1,b-a)}' $d/F.fa; fi
done | awk '{n=length($0); for(i=1;i+24<=n;i+=5){w=substr($0,i,25); if(w!~/[^ACGT]/)print w}}' | sort -u > cand.txt
awk '{print ">"NR; print}' cand.txt > cand.fa
echo "candidates $(wc -l < cand.txt)"
s=$(( $(wc -l < cand.txt) * 2 + 1000000 ))
for k in 1 2; do
	if [ $k = 1 ]; then GN=$T1; else GN=$T2; fi
	jellyfish count -m 25 -C -s $s -t 16 --if cand.fa -o c$k.jf "$GN"
	jellyfish dump -c c$k.jf > c$k.dump
done
# once at most in each haplotype, no low complexity
awk 'FNR==NR{a[$1]=$2; next} $2<=1 && a[$1]<=1{print $1}' c1.dump c2.dump \
	| awk '{w=$0; t=gsub(/A/,"A",w); c=gsub(/C/,"C",w); g=gsub(/G/,"G",w); u=25-t-c-g; if(t>13||c>13||g>13||u>13)next; if($0~/AAAAAA|CCCCCC|GGGGGG|TTTTTT/)next; print}' > words.txt
echo "words $(wc -l < words.txt)"
"$SCF" -s -m 3 -t 16 words.txt "$READS" pool.fa
samtools faidx pool.fa
echo "pool reads $(wc -l < pool.fa.fai)"
rm -f c1.jf c2.jf c1.dump c2.dump cand.fa cand.txt
for d in w/*; do
	. $d/meta
	for x in T F O; do
		if [ -s $d/$x.fa ]; then awk -v p="$id|$x" '/^>/{print ">"p"|"substr($0,2); next}{print}' $d/$x.fa; fi
	done
done > allwin.fa
minimap2 -x map-hifi -N 20 -t 32 allwin.fa pool.fa 2>/dev/null | awk -F'\t' '$4-$3>=500{split($6,p,"|"); print p[1]"\t"$1}' | sort -u > hits.tsv
for d in w/*; do : > $d/names; done
awk -F'\t' '{f="w/"$1"/names"; print $2 >> f; close(f)}' hits.tsv
echo "allotted"
awk -F'\t' '{n[$1]++} END{for(i in n)print n[i]}' hits.tsv | sort -n | awk '{a[NR]=$1}END{print "reads per item min",a[1],"median",a[int((NR+1)/2)],"max",a[NR]}'
