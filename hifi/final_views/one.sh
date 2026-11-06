#!/bin/bash
# one item: phase calls, reads on the view, alignment, BED. Argument: item id
set -uo pipefail
HERE=$(dirname "$(readlink -f "$0")")
. "$HERE/cfg.sh"
cd "$W/w/$1" || exit 1
. ./meta
O=$W/out/$id; mkdir -p "$O"
if [ $a = dep1 ]; then oh=hap2; else oh=hap1; fi
samtools faidx "$W/pool.fa" -r names > reads.fa 2>/dev/null
map() { minimap2 -x map-hifi -c --cs=short --secondary=no -t 2 "$1" reads.fa 2>/dev/null > "$2"; }
map T.fa T.paf; map F.fa F.paf
if [ -s O.fa ]; then map O.fa O.paf; else : > O.paf; fi
rm -f rows.tsv
awk -F'\t' -v OFS='\t' -v g0=$g0 -v g1=$g1 -v f0=$f0 -v f1=$f1 -v id="$id" -v tag=view \
	-v Tf=T.paf -v Of=O.paf -v Ff=F.paf -v rows=rows.tsv -f "$HERE/phase_cmp.awk" T.paf O.paf F.paf > sum.tsv
# rows: id tag read call shared dt dd x cov xs ind bp big span. Label: hap1, hap2, tie, nophase, with _cross for a read crossing the gap
awk -F'\t' -v h=$h -v o=$oh '{c=$4; n=(c=="T")?h:(c=="O")?o:(c=="tie")?"tie":"nophase"; print $3, n ($8==1?"_cross":"")}' rows.tsv > calls.txt
view() {
	minimap2 -ax map-hifi --secondary=no -t 2 "$1" reads.fa 2>/dev/null | samtools view -h -F 0x904 -q 20 -e 'rlen>=1000' - \
		| awk -v calls=calls.txt -f "$HERE/clip.awk" | samtools sort -o "$2" -
	samtools index "$2"
}
view V.fa "$O/view.bam"
cp V.fa "$O/view.fa"
if [ "$VW" = 25000 ] && [ -s O.fa ]; then
	awk -v o=$oh '/^>/{n++; print ">"o (n>1?"_"n:""); next}{print}' O.fa > "$O/other.fa"
	view "$O/other.fa" "$O/other.bam"
fi
# alignment rows
rm -f bet.loc vf.span
{
	cat D.fa
	if [ "$closed" = 1 ]; then cat V.fa; fi
	if [ -s O.fa ]; then
		awk -v n=$v0 '!/^>/{print ">L"; print substr($0,1,n)}' V.fa > vf.fa
		awk -v n=$v1 '!/^>/{print ">R"; print substr($0,n+1)}' V.fa >> vf.fa
		minimap2 -x asm20 --secondary=no -t 2 O.fa vf.fa 2>/dev/null | cut -f1-12 > vf.paf
		awk -F'\t' '{ if(!($1 in M)||$10>M[$1]){M[$1]=$10; T[$1]=$6; A[$1]=$8; B[$1]=$9} }
			END{ if(("L" in T)&&("R" in T)&&T["L"]==T["R"]){ a=A["L"]<A["R"]?A["L"]:A["R"]; b=B["L"]>B["R"]?B["L"]:B["R"]; print T["L"],a,b }
			     else { if("L" in T)print T["L"],A["L"],B["L"]; if("R" in T)print T["R"],A["R"],B["R"] } }' vf.paf > vf.span
		st=$(awk -F'\t' '{m[$5]+=$10} END{print (m["-"]>m["+"])?"-":"+"}' vf.paf)
		# gene items keep the whole other window so that every exon is in the row
		if [ "$VW" = 25000 ]; then awk '/^>/{t=substr($0,2); next}{print t,0,length($0)}' O.fa > vf.span; fi
		n=0
		while read -r t x y; do
			n=$(( n + 1 ))
			awk -v t=">$t" -v x=$x -v y=$y -v st=$st -v o=$oh -v n=$n -f "$HERE/rc.awk" \
				-f <(printf '%s\n' '$0==t{ getline q; s=substr(q,x+1,y-x); print ">"o (n>1?"_"n:""); print (st=="-")?rc(s):s }') O.fa
		done < vf.span
	fi
	b=$(awk -F'\t' -v i=$id '$1==i{print $2"\t"$3}' "$W/bet.tsv")
	if [ -n "$b" ]; then
		reg=${b%%$'\t'*}; st=${b##*$'\t'}
		samtools faidx "$BET" "$reg" | awk -v st=$st -f "$HERE/rc.awk" \
			-f <(printf '%s\n' '!/^>/{s=s toupper($0)} END{print ">bet_primary"; print (st=="-")?rc(s):s}')
		echo "$reg $st" > bet.loc
	fi
	# reads that cross, longest first: up to 3 of this haplotype, 1 of the other
	samtools view "$O/view.bam" | awk -F'\t' -v h=$h -v o=$oh '{ split($1,p,"|"); if(p[1]==h"_cross"||p[1]==o"_cross") print p[1],length($10),$10 }' \
		| sort -k1,1 -k2,2nr \
		| awk -v h=$h -v o=$oh '$1==h"_cross"&&a<3{a++; print ">read_"h"_"a; print $3} $1==o"_cross"&&b<1{b++; print ">read_"o"_"b; print $3}'
} > rows.fa
timeout 3000 mafft --thread 2 --nuc --retree 2 --maxiterate 0 --quiet rows.fa > "$O/aln.fa" 2> mafft.log
# BED, ungapped positions on each row
{
	printf '%s_deposited\t%s\t%s\tN run, %s bp\t0\t+\t%s\t%s\t200,0,0\n' $h $d0 $d1 $(( d1 - d0 )) $d0 $d1
	if [ "$closed" = 1 ]; then
		e1=$v1; if [ $v1 -eq $v0 ]; then e1=$(( v0 + 1 )); fi
		printf '%s\t%s\t%s\tfill, %s bp\t0\t+\t%s\t%s\t200,120,0\n' $vn $v0 $e1 $(( v1 - v0 )) $v0 $e1
	fi
} > base.bed
: > exons.bed
case $kind in *gene*|*control*)
	: > exons.fa
	for gid in ${lab//,/ }; do
		g=${gid%%_*}
		ex=$(awk -F'\t' -v i=$gid '$1==i{print $2}' $GL/prep/exons.tsv)
		awk -v k="$gid|reg|10000" -v g=$g -v ex="$ex" '
			/^>/{p=(substr($0,2)==k); next} p{s=s $0}
			END{ n=split(ex,e,","); s=toupper(s); for(i=1;i<=n;i++){ split(e[i],q,"-"); print ">"g"_exon"i"_"q[2]-q[1]; print substr(s,q[1]+1,q[2]-q[1]) } }' $GL/prep/query.fa >> exons.fa
	done
	awk '/^>read_/{p=0;next} /^>/{p=1} p' rows.fa > rows_noread.fa
	ssearch36 -n -m 8 -E 1 -T 2 exons.fa rows_noread.fa 2>/dev/null | awk -f "$HERE/exon_hits.awk" | sort -k1,1 -k2,2n > exons.bed
	;; esac
cat base.bed exons.bed > "$O/aln.bed"
awk -F'\t' -v n=$vn '$1==n' base.bed exons.bed > "$O/view.bed"
if [ -f "$O/other.fa" ]; then awk -F'\t' -v o=$oh '$1==o||index($1,o"_")==1' exons.bed > "$O/other.bed"; fi
cp sum.tsv rows.tsv meta "$O/"
if [ -f bet.loc ]; then cp bet.loc "$O/"; fi
if [ -f vf.span ]; then cp vf.span "$O/other.span"; fi
echo "$id T_reads/O_reads/T_cross/O_cross/verdict $(cut -f2,3,8,9,21 sum.tsv | tr '\t' ' ') view $(samtools view -c "$O/view.bam") aln $(wc -c < "$O/aln.fa")"
