#!/bin/bash
# Mark every fill in the repaired assemblies with its phase verdict. Sequences are not changed.
# Writes repair/depN.{genesearch,repaired}.phase.bed: chrom start end gap|verdict|n_target_reads|median_substitutions 0 + start end colour
set -euo pipefail
G=/home/toki/gof; F=$G/phase_test/bash_all/fills.tsv
cd "$G/repair" || exit 1
export LC_ALL=C
for a in dep1 dep2; do
	for v in genesearch repaired; do
		awk -F'\t' -v OFS='\t' -v a="$a" '
			BEGIN {
				col["wrong_allele"] = "200,0,0"; col["wrong_allele_structure"] = "200,0,0"
				col["own_or_identical"] = "0,140,0"; col["undecided"] = "220,160,0"
				col["no_target_cross"] = "120,120,120"; col["flanks_tie"] = "120,120,120"; col["no_phasable_reads"] = "120,120,120"
				col["no_fill_sequence"] = "200,200,200"
			}
			FILENAME == ARGV[1] { if ($1 == a) { vd[$2] = $NF; n[$2] = $22; m[$2] = $24 } next }
			{
				split($4, p, "|"); g = p[1]
				w = (g in vd) ? vd[g] : "no_fill_sequence"
				print $1, $2, $3, g "|" w "|" ((g in vd) ? n[g] : "NA") "|" ((g in vd) ? m[g] : "NA"), 0, "+", $2, $3, col[w]
			}' "$F" "$a.$v.bed" > "$a.$v.phase.bed"
		echo "$a.$v.phase.bed $(wc -l < "$a.$v.bed") in, $(wc -l < "$a.$v.phase.bed") out"
		awk -F'\t' '{split($4,p,"|"); print "   "p[2]}' "$a.$v.phase.bed" | sort | uniq -c | sort -rn
	done
done
grep -P '^OZ418347.1\t' dep1.genesearch.phase.bed | grep gap0441
rm -f $G/phase_test/pool.fq.gz $G/phase_test/tests.pkl
du -sh $G/phase_test
