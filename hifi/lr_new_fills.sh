#!/bin/bash
# Fills lr found for gaps still N in the latest assembly (r2_hicanu genesearch2), each put through the phase test.
# Stage 1 list, stage 2 other-haplotype window from the flanks, stage 3 per fill: reads mapped alone to T, O, F, phase_cmp.awk.
set -uo pipefail
G=/home/toki/gof; W=$G/lr_all; P=$G/phase_test/bash_all; R=$G/asm2/r2_hicanu
T1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
T2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
export LC_ALL=C
export PATH=/home/toki/miniforge3/envs/sinederella/bin:/home/toki/miniforge3/envs/hifiasm/bin:$PATH
cd "$W" || exit 1
mkdir -p nf
: > nf/list.tsv
for d in dep1 dep2; do
	if [ $d = dep1 ]; then GT=$T1; GO=$T2; else GT=$T2; GO=$T1; fi
	# still N after round 2 = in repair gap list of the latest genesearch: not in changes2
	awk -F'\t' -v OFS='\t' -v d=$d '
		FILENAME==ARGV[1]{ if(FNR>1) done[$2]=1; next }
		FNR==1{next}
		!($5 in done) && ($9=="insert"||$9=="join"||$9=="genesearch"){ print d,$1,$2,$3,$4,$5,$8,$9,$10,$11,$12,$13 }' \
		$R/$d.changes2.tsv $d.vs_repair.tsv >> nf/list.tsv
done
echo "new fills: $(wc -l < nf/list.tsv)"; cut -f1,8 nf/list.tsv | sort | uniq -c
# flanks of each new fill on the other haplotype
for d in dep1 dep2; do
	if [ $d = dep1 ]; then GT=$T1; GO=$T2; else GT=$T2; GO=$T1; fi
	awk -F'\t' -v d=$d '$1==d{a=$4-15000; if(a<1)a=1; print $2"\t"$3":"a"-"$4-1"\tL"; print $2"\t"$3":"$5+1"-"$5+15000"\tR"}' nf/list.tsv > nf/$d.q
	cut -f2 nf/$d.q > nf/$d.reg
	samtools faidx "$GT" -r nf/$d.reg | awk '/^>/{if(n++)printf "\n"; next}{printf "%s",$0} END{printf "\n"}' > nf/$d.seq
	paste nf/$d.q nf/$d.seq | awk -F'\t' '{print ">"$1"|"$3"\n"toupper($4)}' > nf/$d.flanks.fa
	minimap2 -x asm20 --secondary=no -t 24 "$GO" nf/$d.flanks.fa 2>/dev/null | cut -f1-12 > nf/$d.paf
	# best hit per flank; other window spans both when on one sequence and within 200 kb, else each flank alone
	awk -F'\t' -v OFS='\t' '{ if(!($1 in B)||$10>M[$1]){B[$1]=$6"\t"$8"\t"$9; M[$1]=$10} }
		END{ for(k in B){ split(k,p,"|"); g[p[1]]=1 }
		     for(x in g){ l=B[x"|L"]; r=B[x"|R"]; split(l,a,"\t"); split(r,b,"\t");
		       if(l!=""&&r!=""&&a[1]==b[1]){ lo=a[2]<b[2]?a[2]:b[2]; hi=a[3]>b[3]?a[3]:b[3]; if(hi-lo<=200000){print x,a[1]":"lo+1"-"hi; continue} }
		       if(l!="")print x,a[1]":"a[2]+1"-"a[3]; if(r!="")print x,b[1]":"b[2]+1"-"b[3] } }' nf/$d.paf > nf/$d.owin
done
export G W T1 T2
cut -f1,2 nf/list.tsv | tr '\t' ' ' | xargs -P 8 -L 1 bash /home/toki/gof/lr_all/nf_one.sh
head -1 $P/fills.tsv | cut -f12- > nf/sumhead
{ printf 'assembly\tgap_id\tcontig\tstart1\tend1\trepair_gap\tlr_status\taction\tcause\town\tfill_len\tsupport\t'; cat nf/sumhead; } > new_fills.tsv
awk -F'\t' -v OFS='\t' 'FILENAME==ARGV[1]{K[$1,$2]=$0; next} { k=$1 SUBSEP $2; f="nf/"$1"_"$2".sum"; s=""; if((getline s < f)>0) print K[k],s; close(f) }' nf/list.tsv nf/list.tsv >> new_fills.tsv
echo "== verdicts"; awk -F'\t' 'NR>1{print $1"\t"$8"\t"$NF}' new_fills.tsv | sort | uniq -c
echo finished
