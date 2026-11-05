#!/bin/bash
# phase test of every accepted fill. One xargs job per fill, each with its own files. Stops if TMED4 is not wrong_allele.
set -euo pipefail
E=/home/toki/miniforge3/envs
S=$E/sinederella/bin/samtools
G=/home/toki/gof; P=$G/phase_test; B=$P/bash_all
cd "$B" || exit 1
export LC_ALL=C
date
rm -rf out run names; mkdir -p out run names
# tests.tsv: accepted fills with a sequence, joined to the interval of the insert in the built assembly
awk -F'\t' -v OFS='\t' '
	FILENAME ~ /dep1\.genesearch\.bed$/ { split($4, a, "|"); bed["dep1", a[1]] = $2 "\t" $3; next }
	FILENAME ~ /dep2\.genesearch\.bed$/ { split($4, a, "|"); bed["dep2", a[1]] = $2 "\t" $3; next }
	FNR == 1 { next }
	$14 == "yes" && $8 + 0 > 0 {
		k = $1 SUBSEP $2
		if (!(k in bed)) { print "no bed for", $1, $2 > "/dev/stderr"; exit 1 }
		split($10, c, "|"); split(bed[k], z, "\t")
		if (z[2] - z[1] != $8) { print "bed length differs", $1, $2 > "/dev/stderr"; exit 1 }
		print $1 "_" $2, $1, $2, $3, $4, $5, $6, $7, $8, c[2], c[3], c[4], bed[k]
	}' $G/repair/dep1.genesearch.bed $G/repair/dep2.genesearch.bed $G/repair/repairs.tsv > tests.tsv
wc -l tests.tsv
cat $P/pool_dep1.tsv $P/pool_dep2.tsv | sort -u | awk -F'\t' '{f="names/"$2; print $1 >> f; close(f)}'
ls names | wc -l
n=$(wc -l < names/dep1_gap0441); echo "TMED4 names $n"
if [ "$n" -ne 58 ]; then echo "GUARD FAILED: TMED4 should have 58 reads"; exit 1; fi
zcat $P/pool.fq.gz > pool.fq
$S fqidx pool.fq
date
# guard
bash fill_one.sh dep1_gap0441
v=$(awk -F'\t' '{print $NF}' out/dep1_gap0441.tsv)
echo "TMED4 $(cat out/dep1_gap0441.tsv)"
if [ "$v" != wrong_allele ]; then echo "GUARD FAILED: TMED4 is $v"; exit 1; fi
date
cut -f1 tests.tsv | grep -v '^dep1_gap0441$' | xargs -P 12 -n 1 bash fill_one.sh
date
printf 'assembly\tgap\tchrom\tstart\tend\tcause\taction\tins_len\tsource_name\tsource\tstrand\treads\tT_reads\tO_reads\ttie\tforeign\tshort\tone_side\tT_cross\tO_cross\ttie_cross\tT_span\tT_x\tT_x_median\tT_indel_median\tT_big\tO_x_median\ttie_x_median\tT_enter\tT_enter_x_per_kb\tO_enter_x_per_kb\tverdict\n' > fills.tsv
cat out/*.tsv | sort -t$'\t' -k1,1 -k2,2 >> fills.tsv
cat out/*.rows > fills.reads.tsv
echo "rows $(( $(wc -l < fills.tsv) - 1 )) of $(wc -l < tests.tsv)"
echo "fills with reads=0: $(awk -F'\t' 'NR>1 && $12==0' fills.tsv | wc -l)"
echo "== verdicts"; awk -F'\t' 'NR>1{print $NF}' fills.tsv | sort | uniq -c | sort -rn
# compare with the fills Python finished before it stopped
awk -F'\t' '
	FNR == NR { for (i = 12; i <= NF; i++) b[$1, $2, i] = $i; v[$1, $2] = $NF; next }
	FNR == 1 { next }
	{
		k = $1 SUBSEP $2; n++
		if (!(k in v)) { miss++; next }
		p = $NF; if (p == "no_reads") p = "no_phasable_reads"
		if (p != v[k]) { dv++; if (dv <= 10) print "verdict differs", $1, $2, "python", p, "awk", v[k] }
		for (i = 12; i < NF; i++) if (i != 26 && $i + 0 != b[$1, $2, i] + 0) { dc++; if (dc <= 10) print "column differs", $1, $2, i, "python", $i, "awk", b[$1, $2, i] }
		if ($23 != b[$1, $2, 23]) dx++
	}
	END { print "python fills", n, "not in awk", miss + 0, "| verdicts differ", dv + 0, "| number columns differ", dc + 0, "| substitution lists differ", dx + 0 }
' fills.tsv $P/fills_python_partial.tsv
rm -rf pool.fq pool.fq.fai run
date
echo finished
