#!/bin/bash
# All gaps of both deposited trizona haplotypes through sineclose lr 0.6.0: fresh search of the whole HiFi file, every assembly
# as fill source, k-mer copy number. No patched genome yet: the table decides first. Compared gap by gap with repair/repairs.tsv.
set -uo pipefail
G=/home/toki/gof; W=$G/lr_all; S=/home/toki/sineclose_bash
T1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
T2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
export LC_ALL=C
export PATH=/home/toki/miniforge3/envs/rapidgz/bin:/home/toki/miniforge3/envs/sinederella/bin:/home/toki/miniforge3/envs/asm25/bin:/home/toki/miniforge3/envs/lrtools/bin:/home/toki/miniforge3/envs/jf/bin:$PATH
mkdir -p "$W"; cd "$W" || exit 1
[ -s dep1.gaps.tsv ] || python3 $S/sineclose.py scan $T1 -o dep1.gaps.tsv
[ -s dep2.gaps.tsv ] || python3 $S/sineclose.py scan $T2 -o dep2.gaps.tsv
echo "gaps: dep1 $(( $(wc -l < dep1.gaps.tsv) - 1 )), dep2 $(( $(wc -l < dep2.gaps.tsv) - 1 ))"
ALTS="a25h1=$G/asm25/tri.hic.hap1.p_ctg.fa a25h2=$G/asm25/tri.hic.hap2.p_ctg.fa wgsp=$G/wgs/tri.bp.p_ctg.fa wgsa=$G/wgs/tri.bp.a_ctg.fa hicanu=$G/asm2/hicanu/tri.contigs.fasta"
run() { # name genome other
	bash $S/lr_bash/sineclose_lr.sh "$2" "$1.gaps.tsv" --reads hifi=$G/reads/ERR16751478.fastq.gz --alt other="$3" $ALTS \
		--scfish2 $S/fast3/scfish2 -t 24 --genesearch --meryl $G/sqc/hifi.meryl --calib-windows 2000 -o "$1" > "$1.log" 2>&1
	echo "$1 exit $?"
}
date
run dep1 $T1 $T2 &
run dep2 $T2 $T1 &
wait
date
for d in dep1 dep2; do
	tail -1 $d.log | cut -c1-200
	awk -F'\t' 'NR>1{printf "%s ",$0} END{print ""}' /dev/null
	awk -F'\t' 'NR>1{p=$1; if(q!=""&&0)print} {if(NR>1)printf "   %8.1f s  %s\n",$1-p0,s; p0=$1;s=$2}' $d/times.txt
	# gap by gap against the repairs (repairs.tsv start is 0-based)
	awk -F'\t' -v OFS='\t' -v a=$d '
		FILENAME==ARGV[1]{ if($1==a) R[$3,$4+1]=$7 "|" $14 "|" $2; next }
		FNR==1{ print "gap_id","contig","start1","end1","repair_gap","repair_action","repair_accepted","lr_status","lr_action","lr_cause","lr_own","lr_fill_len","lr_support","supported"; next }
		{ k=$2 SUBSEP $3; split((k in R)?R[k]:"none|none|none",r,"|"); print $1,$2,$3,$4,r[3],r[1],r[2],$5,$6,$7,$8,$9,$11,$28 }' \
		$G/repair/repairs.tsv $d/gaps_class.tsv > $d.vs_repair.tsv
	echo "== $d: repair state against lr state"
	awk -F'\t' 'NR>1{ rs=($7=="yes")?(($6=="genesearch")?"rep:genesearch":"rep:strict"):"rep:N"; ls=($9=="insert"||$9=="join")?"lr:strict":($9=="genesearch")?"lr:genesearch":"lr:N"; c[rs" "ls]++ } END{for(k in c)print c[k],k}' $d.vs_repair.tsv | sort -k2,2 -k3,3
done
echo finished
