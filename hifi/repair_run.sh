#!/bin/bash
# Gap repair of both deposited haplotypes (repair_gaps.py): dep2 gap causes -> plan -> unique junction words ->
# one scfish pass over the HiFi reads -> verify each junction -> repaired / gene-search FASTAs.
# Usage: repair_run.sh OUT_DIR [toy]   (toy: dep1 only, gaps in $TOY_GAPS)
set -euo pipefail
O=$1
MODE=${2:-full}
mkdir -p "$O"
cd "$O" || exit 1
C=/home/toki/gof/code/hifi
G=/home/toki/gof/gf_all
R=/home/toki/gof/reads/ERR16751478.fastq.gz
D1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
D2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
TOY_GAPS=${TOY_GAPS:-gap0300,gap0549,gap0551,gap0011,gap0001}
export PATH=/home/toki/miniforge3/envs/hifiasm/bin:/home/toki/miniforge3/envs/jf/bin:/home/toki/miniforge3/envs/sinederella/bin:/home/toki/miniforge3/envs/asm25/bin:/home/toki/miniforge3/envs/lrtools/bin:$PATH
T=32

if [ "$MODE" = toy ]; then
  ASMS="dep1"
  EXTRA="--gaps $TOY_GAPS"
else
  ASMS="dep1 dep2"
  EXTRA=""
  echo "#### dep2 gap causes $(date '+%F %T')"
  mkdir -p causes2
  [ -s "$D2.fai" ] || samtools faidx "$D2"
  python3 $C/gap_flank_cn.py prep "$D2" $G/dep2/gaps.tsv causes2/flanks.fa
  meryl-lookup -wig-count -sequence causes2/flanks.fa -mers /home/toki/gof/sqc/hifi.meryl -output causes2/flanks.wig > causes2/meryl.log 2>&1
  python3 $C/gap_flank_cn.py table causes2/flanks.fa causes2/flanks.wig /home/toki/gof/ends/kmer/gc_calib.txt $G/dep2/gaps.tsv $G/fill_cn.tsv dep1 > causes2/gap_flank_cn.tsv
  (cd causes2 && trf flanks.fa 2 7 7 80 10 50 500 -h -ngs > flanks.trf.dat || true)
  python3 $C/gap_classify.py $G/dep2/gaps.tsv causes2/gap_flank_cn.tsv causes2/flanks.fa causes2/flanks.trf.dat no_X_in_hap2 > causes2/gap_causes.tsv 2> causes2/gap_causes.counts
  rm -f causes2/flanks.wig
fi

for a in $ASMS; do
  if [ $a = dep1 ]; then GF=$D1; OTHER=dep2; CA=$G/flankcn/gap_causes.tsv; else GF=$D2; OTHER=dep1; CA=causes2/gap_causes.tsv; fi
  echo "#### plan $a $(date '+%F %T')"
  python3 $C/repair_gaps.py plan $GF $G/$a/gaps.tsv $G/$a/fills.fa $G/all_fills.fa $G/fill_cn.tsv $CA $OTHER $a $EXTRA
  echo "#### words $a $(date '+%F %T')"
  n=$(grep -c '>' $a/cand.fa)
  jellyfish count -m 25 -C -s $((4 * n + 1000)) -t 16 --if $a/cand.fa -o $a/cand.jf $GF $a/pieces.fa
  jellyfish dump -c $a/cand.jf > $a/cand.counts
  python3 $C/repair_gaps.py words $a $a/cand.counts
  rm -f $a/cand.jf $a/cand.counts
done

echo "#### fish $(date '+%F %T')"
cat $(for a in $ASMS; do echo $a/anchors.txt; done) | sort -u > anchors.txt
wc -l anchors.txt
pigz -dc -p 8 $R | paste - - - - | /home/toki/gof/scfish/scfish anchors.txt 3 \
  | awk -F'\t' '{split($1,h," "); print "@" substr(h[1],2) "\n" $2 "\n+\n" $4}' > fished.fq
echo "fished reads: $(( $(wc -l < fished.fq) / 4 ))"

for a in $ASMS; do
  if [ $a = dep1 ]; then GF=$D1; else GF=$D2; fi
  echo "#### verify $a $(date '+%F %T')"
  python3 $C/repair_gaps.py verify $GF $a fished.fq $a $T $([ "$MODE" = toy ] && echo --edited-only)
done
echo "#### done $(date '+%F %T')"
for a in $ASMS; do python3 $C/repair_gaps.py bed $a $a; done
