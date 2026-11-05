#!/bin/bash
# gene status on the final assemblies, same 191 genes and code as geneloss2; orig1/orig2 calls reused. Waits for final_build.
set -uo pipefail
G=/home/toki/gof
while ! grep -q '^finished' $G/lr_all/final_build.log; do
	ps -u toki -o args | grep -v grep | grep -q final_build.sh || break
	sleep 15
done
grep -q '^finished' $G/lr_all/final_build.log || { echo "final build did not finish"; exit 1; }
mkdir -p $G/geneloss3/prep; cd $G/geneloss3 || exit 1
export PATH=/home/toki/miniforge3/envs/hifiasm/bin:/home/toki/toga2/env_mlc2/bin:/home/toki/miniforge3/envs/sinederella/bin:$PATH
export TMPDIR=/home/toki/tmp
C=$G/code; A=$G/gf_all; R=$G/repair; F=$G/final
D1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
D2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
O=prep
cp $G/geneloss2/prep/genes.tsv $G/geneloss2/prep/query.fa $G/geneloss2/prep/query.bed $G/geneloss2/prep/exons.tsv $G/geneloss2/prep/orig1.calls.tsv $G/geneloss2/prep/orig2.calls.tsv $O/
ls $G/geneloss2/prep | grep -v -E 'paf$|calls|mm.log|final|summary' | head
for a in dep1 dep2 dep1gs dep2gs; do
	case $a in
		dep1) X=$F/dep1.strict.fa; OR=$D1; GP=$A/dep1/gaps.tsv; CA=$A/flankcn/gap_causes.tsv;;
		dep2) X=$F/dep2.strict.fa; OR=$D2; GP=$A/dep2/gaps.tsv; CA=$R/causes2/gap_causes.tsv;;
		dep1gs) X=$F/dep1.genesearch.fa; OR=$D1; GP=$A/dep1/gaps.tsv; CA=$A/flankcn/gap_causes.tsv;;
		dep2gs) X=$F/dep2.genesearch.fa; OR=$D2; GP=$A/dep2/gaps.tsv; CA=$R/causes2/gap_causes.tsv;;
	esac
	echo "#### $a $(date '+%T')"
	minimap2 -c --cs -x asm20 -N 5 -t 32 $X $O/query.fa > $O/$a.paf 2> $O/$a.mm.log
	python3 $C/hifi/geneloss_repaired.py call $O $a $O/$a.paf $X $OR $GP $CA || { echo "call failed $a"; exit 1; }
done
python3 $C/hifi/geneloss_repaired.py combine $O
cp $O/final.tsv $O/summary.txt .
rm -f $O/*.paf
echo "#### old against new, per haplotype call"
for k in 1 2; do :; done
head -1 final.tsv | tr '\t' '\n' | cat -n | head -16
cut -c1-200 summary.txt | head -20
echo finished
