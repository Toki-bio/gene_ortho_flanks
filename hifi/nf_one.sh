#!/bin/bash
# phase test of one lr fill. Args: assembly gap_id. Reads = the reads lr fished for this haplotype, mapped alone to each window.
set -uo pipefail
a=$1; g=$2; FL=25000
if [ "$a" = dep1 ]; then GT=$T1; GO=$T2; else GT=$T2; GO=$T1; fi
L=$W/lr_all_lib
cd "$W/nf" || exit 1
d=w_${a}_$g; mkdir -p $d
IFS=$'\t' read -r _ _ c s e _ _ _ _ own _ _ < <(awk -F'\t' -v a=$a -v g=$g '$1==a&&$2==g' list.tsv)
IFS=$'\t' read -r _ _ ty val _ < <(awk -F'\t' -v o="$own" -v g=$g '$1==g&&$2==o' ../$a/cands.tsv)
l=$(( s - 1 < FL ? s - 1 : FL ))
samtools faidx "$GT" "$c:$(( s - l ))-$(( s - 1 ))" | grep -v '>' | tr -d '\n' | tr 'acgtn' 'ACGTN' > $d/L
samtools faidx "$GT" "$c:$s-$e" | grep -v '>' | tr -d '\n' | tr 'acgtn' 'ACGTN' > $d/H
samtools faidx "$GT" "$c:$(( e + 1 ))-$(( e + FL ))" | grep -v '>' | tr -d '\n' | tr 'acgtn' 'ACGTN' > $d/R
{ echo ">T"; cat $d/L $d/H $d/R; echo; } > $d/T.fa
ov=0; fill=""
if [ "$ty" = I ]; then ov=$(( -val )); else fill=$val; fi
{ echo ">F"; cat $d/L; printf '%s' "$fill"; tail -c +$(( ov + 1 )) $d/R; echo; } > $d/F.fa
fl=${#fill}
awk -F'\t' -v g=$g '$1==g' FS=' ' OFS='\t' /dev/null
: > $d/oreg
awk -v g=$g '$1==g{print $2}' FS='\t' $a.owin > $d/oreg
[ -s $d/oreg ] || { printf '0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t\tNA\tNA\t\tNA\tNA\t0\tNA\tNA\tno_other_window\n' > ${a}_$g.sum; rm -rf $d; exit 0; }
samtools faidx "$GO" -r $d/oreg | awk '/^>/{n++; print ">O"n; next}{print}' > $d/O.fa
FQ=$W/$a/fished_hifi.fa
map() { minimap2 -x map-hifi -c --cs=short --secondary=no -t 2 "$1" "$FQ" 2>/dev/null > "$2"; }
map $d/T.fa $d/T.paf; map $d/O.fa $d/O.paf; map $d/F.fa $d/F.paf
h=$(wc -c < $d/H)
awk -F'\t' -v OFS='\t' -v g0=$l -v g1=$(( l + h )) -v f0=$l -v f1=$(( l + fl )) -v id="${a}_$g" -v tag=lr \
	-v Tf=$d/T.paf -v Of=$d/O.paf -v Ff=$d/F.paf -v rows=$d/rows -f $L/lib.awk -f $L/phase_cmp_lr.awk $d/T.paf $d/O.paf $d/F.paf > ${a}_$g.sum
rm -rf $d
