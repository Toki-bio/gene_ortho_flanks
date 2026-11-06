#!/bin/bash
# windows per item: T deposited with N, F final with fill, V view cut, D deposited view cut, O other haplotype from the flanks
set -uo pipefail
. "$(dirname "$0")/cfg.sh"
cd "$W" || exit 1
L=${1:-items.tsv}
up() { grep -v '>' | tr -d '\n' | tr 'acgtn' 'ACGTN'; }
mkdir -p w
: > fl.dep1.fa; : > fl.dep2.fa
while IFS=$'\t' read -r id kind a gap c s e lab st rest; do
	d=w/$id; mkdir -p $d; rm -f $d/oreg
	if [ $a = dep1 ]; then GT=$T1; h=hap1; else GT=$T2; h=hap2; fi
	len=$(awk -F'\t' -v c=$c '$1==c{print $2}' $GT.fai)
	ta=$(( s - FL < 1 ? 1 : s - FL )); tb=$(( e + FL > len ? len : e + FL ))
	{ echo ">T"; samtools faidx $GT $c:$ta-$tb | up; echo; } > $d/T.fa
	g0=$(( s - ta )); g1=$(( e - ta + 1 ))
	VW=5000; case $kind in *gene*|*control*) VW=25000;; esac
	line=$(awk -F'\t' -v g=$gap '{split($4,p,"|")} p[1]==g{print $2,$3}' $F/$a.genesearch.marks.bed)
	if [ -n "$line" ]; then
		read -r m0 m1 <<< "$line"
		flen=$(awk -F'\t' -v c=$c '$1==c{print $2}' $F/$a.genesearch.fa.fai)
		fa=$(( m0 + 1 - FL < 1 ? 1 : m0 + 1 - FL )); fb=$(( m1 + FL > flen ? flen : m1 + FL ))
		{ echo ">F"; samtools faidx $F/$a.genesearch.fa $c:$fa-$fb | up; echo; } > $d/F.fa
		f0=$(( m0 - fa + 1 )); f1=$(( m1 - fa + 1 )); closed=1; vn=${h}_final; src=$d/F.fa
	else
		{ echo ">F"; grep -v '>' $d/T.fa; } > $d/F.fa
		f0=$g0; f1=$g1; closed=0; vn=${h}_deposited; src=$d/T.fa
	fi
	# view window, offsets inside it
	va=$(( f0 - VW < 0 ? 0 : f0 - VW ))
	awk -v a=$va -v b=$(( f1 + VW )) -v n=$vn '!/^>/{print ">"n; print substr($0,a+1,b-a)}' $src > $d/V.fa
	ua=$(( g0 - VW < 0 ? 0 : g0 - VW ))
	awk -v a=$ua -v b=$(( g1 + VW )) -v n=${h}_deposited '!/^>/{print ">"n; print substr($0,a+1,b-a)}' $d/T.fa > $d/D.fa
	printf 'id=%s\nkind=%s\na=%s\nh=%s\ngap=%s\nc=%s\ns=%s\ne=%s\nlab=%s\ng0=%s\ng1=%s\nf0=%s\nf1=%s\nv0=%s\nv1=%s\nd0=%s\nd1=%s\nvn=%s\nVW=%s\nta=%s\nclosed=%s\n' \
		$id "$kind" $a $h $gap $c $s $e "$lab" $g0 $g1 $f0 $f1 $(( f0 - va )) $(( f1 - va )) $(( g0 - ua )) $(( g1 - ua )) $vn $VW $ta $closed > $d/meta
	# 15 kb flanks for the other haplotype
	awk -v id=$id -v a=$g0 -v b=$g1 '!/^>/{ l=a<15000?a:15000; print ">"id"|L"; print substr($0,a-l+1,l); print ">"id"|R"; print substr($0,b+1,15000) }' $d/T.fa >> fl.$a.fa
done < "$L"
for a in dep1 dep2; do
	[ -s fl.$a.fa ] || continue
	if [ $a = dep1 ]; then GO=$T2; else GO=$T1; fi
	minimap2 -x asm20 --secondary=no -t 32 "$GO" fl.$a.fa 2>/dev/null | cut -f1-12 > fl.$a.paf
	awk -F'\t' -v OFS='\t' '{ if(!($1 in B)||$10>M[$1]){B[$1]=$6"\t"$8"\t"$9; M[$1]=$10} }
		END{ for(k in B){ split(k,p,"|"); g[p[1]]=1 }
		     for(x in g){ l=B[x"|L"]; r=B[x"|R"]; split(l,u,"\t"); split(r,v,"\t");
		       if(l!=""&&r!=""&&u[1]==v[1]){ lo=u[2]<v[2]?u[2]:v[2]; hi=u[3]>v[3]?u[3]:v[3]; if(hi-lo<=200000){print x,u[1]":"lo+1"-"hi; continue} }
		       if(l!="")print x,u[1]":"u[2]+1"-"u[3]; if(r!="")print x,v[1]":"v[2]+1"-"v[3] } }' fl.$a.paf > fl.$a.owin
	while IFS=$'\t' read -r id reg; do
		printf '%s\n' "$reg" >> w/$id/oreg
	done < fl.$a.owin
done
for d in w/*; do
	[ -s $d/oreg ] || continue
	. $d/meta; if [ $a = dep1 ]; then GO=$T2; else GO=$T1; fi
	sort -u $d/oreg -o $d/oreg
	# gene items: 30 kb more each side, exons of the other haplotype can lie past the flanks
	if [ "$VW" = 25000 ]; then
		awk -v fai="$GO.fai" 'BEGIN{while((getline l<fai)>0){split(l,p,"	"); N[p[1]]=p[2]}} {split($0,q,":"); split(q[2],r,"-"); a=r[1]-30000; b=r[2]+30000; if(a<1)a=1; if(b>N[q[1]])b=N[q[1]]; print q[1]":"a"-"b}' $d/oreg > $d/oreg2
		mv $d/oreg2 $d/oreg
	fi
	samtools faidx "$GO" -r $d/oreg | awk '/^>/{if(n++)printf "\n"; printf ">O%d\n",n; next}{printf "%s",toupper($0)} END{printf "\n"}' > $d/O.fa
done
# betulina primary for each view, one batch
for d in w/*; do . $d/meta; awk -v id=$id '/^>/{print ">"id; next}{gsub(/N/,""); print}' $d/V.fa; done > vq.fa
minimap2 -x asm20 -c --secondary=no -t 32 "$BET" vq.fa 2>/dev/null | cut -f1-11 > vq.paf
awk -F'\t' -v OFS='\t' '{k=$1 SUBSEP $6 SUBSEP $5; m[k]+=$10; if(!(k in lo)||$8<lo[k])lo[k]=$8; if($9>hi[k])hi[k]=$9}
	END{for(k in m){split(k,p,SUBSEP); if(m[k]>b[p[1]]){b[p[1]]=m[k]; r[p[1]]=p[2]":"lo[k]+1"-"hi[k]"\t"p[3]"\t"m[k]}} for(i in r)print i,r[i]}' vq.paf > bet.tsv
echo "windows done: $(ls w | wc -l), other window: $(ls w/*/O.fa | wc -l), betulina: $(wc -l < bet.tsv)"
