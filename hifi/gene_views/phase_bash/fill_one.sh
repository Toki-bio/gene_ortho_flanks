#!/bin/bash
# phase test of one accepted fill; argument = fill id (dep1_gap0441). Reads tests.tsv, win_dep*.tsv, names/, pool.fq
set -euo pipefail
E=/home/toki/miniforge3/envs
S=$E/sinederella/bin/samtools; M=$E/hifiasm/bin/minimap2
G=/home/toki/gof; P=$G/phase_test; B=$P/bash_all; FL=25000
id=$1
T1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
T2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
export LC_ALL=C
IFS=$'\t' read -r _ asm gap chrom g0 g1 cause action ins src_name src strand s e < <(awk -F'\t' -v i="$id" '$1==i' "$B/tests.tsv")
if [ "$asm" = dep1 ]; then GT=$T1; GO=$T2; other=dep2; else GT=$T2; GO=$T1; other=dep1; fi
d=$B/run/$id
mkdir -p "$d"
cd "$d" || exit 1
l=$(( g0 < FL ? g0 : FL )); lf=$(( s < FL ? s : FL ))
$S faidx "$GT" "$chrom:$(( g0 - l + 1 ))-$(( g1 + FL ))" > T.fa
$S faidx "$G/repair/$asm.genesearch.fa" "$chrom:$(( s - lf + 1 ))-$(( e + FL ))" > F.fa
awk -F'\t' -v i="$id" '$4==i{print $1":"$2+1"-"$3}' "$P/win_$other.tsv" > oreg
[ -s oreg ] || { echo "$id no other window" >&2; exit 1; }
$S faidx "$GO" -r oreg > O.fa
$S fqidx "$B/pool.fq" -r "$B/names/$id" > reads.fq 2> fqidx.log
[ -s reads.fq ] || { echo "$id no reads pulled" >&2; exit 1; }
map() { $M -x map-hifi -c --cs=short --secondary=no -t 2 "$1" reads.fq 2>/dev/null > "$2"; }
map T.fa T.paf; map O.fa O.paf; map F.fa F.paf
h1=$(( l + g1 - g0 )); f1=$(( lf + e - s ))
sum=$(awk -F'\t' -v OFS='\t' -v g0=$l -v g1=$h1 -v f0=$lf -v f1=$f1 -v id="$id" -v tag=built \
	-v Tf=T.paf -v Of=O.paf -v Ff=F.paf -v rows="$B/out/$id.rows" -f "$B/phase_cmp.awk" T.paf O.paf F.paf)
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$asm" "$gap" "$chrom" "$g0" "$g1" "$cause" "$action" "$ins" "$src_name" "$src" "$strand" "$sum" > "$B/out/$id.tsv"
cd "$B" || exit 1
rm -rf "$d"
