#!/bin/bash
# one control site: reads named in the Python table, mapped alone to T, O and each fixed window, then phase_cmp.awk
set -euo pipefail
E=/home/toki/miniforge3/envs
S=$E/sinederella/bin/samtools; M=$E/hifiasm/bin/minimap2
P=/home/toki/gof/phase_test
B=$P/bash_ctl
id=$1
d=$B/$id
mkdir -p "$d"
cd "$d" || exit 1
export LC_ALL=C
awk -F'\t' -v i="$id" '$1==i{print $2}' $P/controls.reads.tsv | sort -u > names
$S fqidx $B/ctl.fq -r names > reads.fq 2> fqidx.log
map() { $M -x map-hifi -c --cs=short --secondary=no -t 2 "$1" reads.fq 2>/dev/null > "$2"; }
map $P/$id/T.fa T.paf
map $P/$id/O.fa O.paf
seqlen() { grep -v '>' "$1" | tr -d '\n' | wc -c; }
tl=$(seqlen $P/$id/T.fa)
for f in $P/$id/F_*.fa; do
	tag=$(basename "$f" .fa); tag=${tag#F_}
	map "$f" "F_$tag.paf"
	fl=$(( $(seqlen "$f") - tl + 100 ))
	f1=$(( 25000 + fl ))
	awk -F'\t' -v OFS='\t' -v g0=25000 -v g1=25100 -v f0=25000 -v f1=$f1 -v id="$id" -v tag="$tag" \
		-v Tf=T.paf -v Of=O.paf -v Ff="F_$tag.paf" -v rows="rows_$tag.tsv" \
		-f $B/phase_cmp.awk T.paf O.paf "F_$tag.paf" > "sum_$tag.tsv"
done
rm -f reads.fq T.paf O.paf F_*.paf
echo "$id done"
