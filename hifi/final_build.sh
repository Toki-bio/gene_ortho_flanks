#!/bin/bash
# Final trizona assemblies: round-2 repaired assemblies plus the new lr fills, four FASTA + marks. Gap ids are matched by sequence,
# because coordinates differ between deposited and repaired. Strict = joins and inserts not proven wrong allele; gene-search = all.
set -uo pipefail
G=/home/toki/gof; W=$G/lr_all; R=$G/asm2/r2_hicanu; O=$G/final; L=/home/toki/sineclose_bash/lr_bash; S=/home/toki/sineclose_bash
T1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
T2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
export LC_ALL=C
export PATH=/home/toki/miniforge3/envs/sinederella/bin:/home/toki/miniforge3/envs/asm25/bin:$PATH
mkdir -p "$O"; cd "$O" || exit 1
one() { # dep version basefa basebed
	d=$1; v=$2; B=$3; BB=$4; t=$d.$v
	GD=$T1; [ $d = dep2 ] && GD=$T2
	[ -s "$B.fai" ] || samtools faidx "$B"
	python3 $S/sineclose.py scan "$B" -o $t.runs.tsv > /dev/null
	tail -n +2 $t.runs.tsv | cut -f1-4 > $t.runs4
	# flank keys of N runs in the base, and of the new fills on the deposited genome
	awk -F'\t' '{print $2":"$3-40"-"$3-1; print $2":"$4+1"-"$4+40}' $t.runs4 > $t.rk
	samtools faidx "$B" -r $t.rk | awk '/^>/{if(n++)printf "\n"; next}{printf "%s",toupper($0)} END{printf "\n"}' | paste - - > $t.rkeys
	paste $t.runs4 $t.rkeys > $t.runs5
	cp $t.runs5 $t.runs5k
	awk -F'\t' -v d=$d -v v=$v 'NR>1 && $1==d { w=$NF; a=$8; if(v=="strict" && !(a=="insert"||a=="join") ) next; if(v=="strict" && w ~ /^wrong_allele/) next; print }' $W/new_fills.tsv > $t.take
	awk -F'\t' '{print $3":"$4-40"-"$4-1; print $3":"$5+1"-"$5+40}' $t.take > $t.dk
	samtools faidx "$GD" -r $t.dk | awk '/^>/{if(n++)printf "\n"; next}{printf "%s",toupper($0)} END{printf "\n"}' | paste - - > $t.dkeys
	paste $t.take $t.dkeys > $t.take2
	# base run id for each new fill; fills and candidate rows under that id
	awk -F'\t' -v OFS='\t' -v f=$t.fills -v m=$t.map '
		FILENAME==ARGV[1]{ B[$2,$5,$6]=$1; next }
		{ n=NF; k=$3 SUBSEP $(n-1) SUBSEP $n; if(k in B){ print B[k],$10,$8 > f; print $2,B[k],$6 > m } else print "NOT MATCHED",$2,$6 > "/dev/stderr" }' $t.runs5 $t.take2
	awk -F'\t' -v OFS='\t' 'FILENAME==ARGV[1]{split($0,x,"\t"); M[x[1]]=x[2]; next}{ if($1 in M){ g=$1; $1=M[g]; print } }' $t.map $W/$d/cands.tsv > $t.cands
	awk -f $L/lib.awk -v plan=$t.plan -v changes=$t.changes.tsv -v bed=$t.newbed -v regions=$t.regions -v shifts=$t.shifts -f $L/patch_plan.awk "$B.fai" $t.runs4 $t.fills $t.cands
	samtools faidx "$B" -r $t.regions > $t.seqs
	awk -v seqf=$t.seqs -v keep=1 -f $L/patch_write.awk $t.plan > $t.fa
	samtools faidx $t.fa
	# marks: base marks moved to final coordinates, plus new fills
	awk -F'\t' -v OFS='\t' -v d=$d -v vf=$G/phase_test/bash_all/fills.tsv -v nf=$W/new_fills.tsv -v mp=$t.map -v sf=$t.shifts '
		function fin(c,p,   i,s){ s=0; for(i=1;i<=n[c];i++){ if(OS[c,i]-1<p) s+=SH[c,i]; else break } return p+s }
		BEGIN{ while((getline l<vf)>0){ split(l,a,"\t"); if(a[1]==d) V[a[2]]=a[32] }
		       while((getline l<mp)>0){ split(l,a,"\t"); ID[a[2]]=a[1]; RG[a[2]]=a[3] }
		       while((getline l<nf)>0){ m=split(l,a,"\t"); if(a[1]==d) NV[a[2]]=a[m] }
		       while((getline l<sf)>0){ split(l,a,"\t"); SF[a[1]]=a[2] } }
		FILENAME==ARGV[1]{ if(FNR>1){ c=$2; i=++n[c]; OS[c,i]=$3; SH[c,i]=SF[$1]; NS[$1]=$7; NE[$1]=$8; AC[$1]=$5; CC[$1]=$2 } next }
		{ split($4,p,"|"); g=p[1]; print $1,fin($1,$2),fin($1,$3),g"|" p[2] "|" p[4] "|" ((g in V)?V[g]:"not_tested") }
		END{ for(x in NS){ g=ID[x]; print CC[x],NS[x]-1,NE[x],RG[x]"|lr|"AC[x]"|"NV[g] } }' $t.changes.tsv "$BB" | sort -k1,1 -k2,2n > $t.marks.bed
	# oracle: every new fill read back from the final FASTA between its flanks; left 40 of the deposited gap, the fill, then the base right flank
	awk -F'\t' -v OFS='\t' 'FILENAME==ARGV[1]{ if(FNR>1){ S[$1]=$7; E[$1]=$8; C[$1]=$2 } next } { x=$1; a=S[x]-41; b=E[x]; print x,C[x]":"a+1"-"S[x]-1,C[x]":"S[x]"-"b }' $t.changes.tsv $t.fills > $t.chk
	cut -f2 $t.chk > $t.ck1
	samtools faidx $t.fa -r $t.ck1 | awk '/^>/{if(n++)printf "\n"; next}{printf "%s",toupper($0)} END{printf "\n"}' > $t.ck1s
	awk -F'\t' 'FILENAME==ARGV[1]{ K[$1]=$(NF-1); next } FILENAME==ARGV[2]{ i++; L[i]=$0; next } { j++; if(L[j]==K[$1])ok++; else bad++ } END{ print "   left flank before each new fill: ok " ok+0 ", wrong " bad+0 }' $t.runs5k $t.ck1s $t.chk
	echo "$t: base runs $(wc -l < $t.runs4), new fills $(wc -l < $t.fills), final runs $(python3 $S/sineclose.py scan $t.fa -o /dev/stdout 2>/dev/null | grep -c '^gap')"
	echo "   base bp $(awk '{s+=$2} END{print s}' "$B.fai"), final bp $(awk '{s+=$2} END{print s}' $t.fa.fai), sum of changes $(awk -F'\t' '{s+=$2} END{print s+0}' $t.shifts)"
	rm -f $t.rk $t.rkeys $t.dk $t.dkeys $t.take2 $t.seqs $t.plan $t.regions $t.newbed $t.chk $t.ck1 $t.ck1s $t.runs5 $t.runs5k $t.runs4
}
for d in dep1 dep2; do
	one $d strict $R/$d.repaired2.fa $R/$d.repaired2.bed
	one $d genesearch $R/$d.genesearch2.fa $R/$d.genesearch2.bed
done
ls -la *.fa | cut -c25-
cut -f4 dep1.genesearch.marks.bed | cut -d'|' -f2,4 | sort | uniq -c | sort -rn | head -14
echo finished
