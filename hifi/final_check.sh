#!/bin/bash
# after final_build and geneloss3: true N-run counts, fills holding N, gene changes, status of every gap
set -uo pipefail
G=/home/toki/gof; F=$G/final; W=$G/lr_all
export LC_ALL=C
export PATH=/home/toki/miniforge3/envs/sinederella/bin:$PATH
cd "$F" || exit 1
echo "== N runs >= 10 and N bases"
for f in dep1.strict dep1.genesearch dep2.strict dep2.genesearch; do
	python3 /home/toki/sineclose_bash/sineclose.py scan $f.fa -o $f.final_runs.tsv > /dev/null
	echo "$f runs $(( $(wc -l < $f.final_runs.tsv) - 1 )), N bp $(awk -F'\t' 'NR>1{s+=$5} END{print s+0}' $f.final_runs.tsv)"
done
echo "== new fills that hold N"
for d in dep1 dep2; do
	awk -F'\t' -v d=$d 'NR==FNR{ if($1==d) O[$2]=$10; next } ($1 in O) && O[$1]==$2 && $3=="S" && toupper($4) ~ /N/ { s=toupper($4); n=gsub(/N/,"",s); print d,$1,$2,length($4),"N bases",n }' $W/new_fills.tsv $W/$d/cands.tsv
done
echo "== genes: old and new per-haplotype calls"
awk -F'\t' -v OFS='\t' 'FILENAME==ARGV[1]{ a[$1]=$6; b[$1]=$15; c[$1]=$7; e[$1]=$16; x[$1]=$10; next }
	FNR>1 { if(a[$1]!=$6) print "hap1 strict",$1,a[$1],"->",$6,x[$1],"->",$10; if(b[$1]!=$15) print "hap1 gs",$1,b[$1],"->",$15; if(c[$1]!=$7) print "hap2 strict",$1,c[$1],"->",$7; if(e[$1]!=$16) print "hap2 gs",$1,e[$1],"->",$16 }' $G/geneloss2/final.tsv $G/geneloss3/final.tsv | sort | head -40
echo "hap1 gene-search still gap:"; awk -F'\t' 'NR>1 && $15 ~ /gap/{printf "%s ",$1} END{print ""}' $G/geneloss3/final.tsv
echo "hap1 strict counts among 150:"; awk -F'\t' 'NR>1 && $3=="assembly_gap"{print $6}' $G/geneloss3/final.tsv | sort | uniq -c | tr '\n' ' '; echo
echo "hap1 gs counts among 150:"; awk -F'\t' 'NR>1 && $3=="assembly_gap"{print $15}' $G/geneloss3/final.tsv | sort | uniq -c | tr '\n' ' '; echo
echo "hap2 gs counts among 150:"; awk -F'\t' 'NR>1 && $3=="assembly_gap"{print $16}' $G/geneloss3/final.tsv | sort | uniq -c | tr '\n' ' '; echo
# every gap of both haplotypes, one row: where it stands in the final assemblies and why
{
printf 'assembly\tgap\tcontig\tstart1\tend1\tfinal_state\tround\taction\tverdict\tlr_status\tlr_cause\n'
for d in dep1 dep2; do
	awk -F'\t' -v OFS='\t' -v d=$d '
		FILENAME==ARGV[1]{ split($4,p,"|"); S[p[1]]=1; next }
		FILENAME==ARGV[2]{ split($4,p,"|"); R[p[1]]=p[2]; A[p[1]]=p[3]; V[p[1]]=p[4]; next }
		FNR==1{next}
		{ g=$5; st=(g in S)?"closed_strict":(g in R)?"closed_genesearch":"N"; print d,g,$2,$3,$4,st,(g in R)?R[g]:"",(g in R)?A[g]:"",(g in R)?V[g]:"",$8,$10 }' \
		$d.strict.marks.bed $d.genesearch.marks.bed $W/$d.vs_repair.tsv
done
} > gap_status.tsv
echo "== final state"; awk -F'\t' 'NR>1{print $1"\t"$6}' gap_status.tsv | sort | uniq -c
echo "== still N, by lr result"; awk -F'\t' 'NR>1 && $6=="N"{print $1"\t"$10"\t"$11}' gap_status.tsv | sort | uniq -c | sort -k2,2 -k1,1nr
echo "== rows without repair gap id: $(awk -F'\t' 'NR>1 && $2==""' gap_status.tsv | wc -l)"
echo finished
