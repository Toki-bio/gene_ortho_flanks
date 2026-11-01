#!/bin/bash
# Repair round 2: a new genome-wide assembly (+ the unpurged hifiasm 0.16.1 -l0 p_ctg/a_ctg) as extra fill sources for
# the gaps that round 1 (repair_run.sh, /home/toki/gof/repair) left in both deposited haplotypes. Same logic as round 1:
# wgs/gapfill.py spans -> fill_cn.py copy number -> repair_gaps.py plan/words -> one scfish pass -> verify (>= 3 crossing
# reads each junction, identity >= 0.97) -> asm2_merge.py (round-1 + round-2 edits on the deposited sequence).
# Usage: asm2_repair2.sh OUT_DIR NEW_ASM.fa|none [toy]   (toy: dep1 only, gaps restricted to TOY_GAPS, reads = round-1 fished)
set -euo pipefail
O=$1
NEW=$2
MODE=${3:-full}
mkdir -p "$O"
cd "$O" || exit 1
C=/home/toki/gof/code/hifi
W=/home/toki/gof/code/wgs
R=/home/toki/gof/reads/ERR16751478.fastq.gz
R1=/home/toki/gof/repair
D1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
D2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
WP=/home/toki/gof/wgs/tri.bp.p_ctg.fa
WA=/home/toki/gof/wgs/tri.bp.a_ctg.fa
export PATH=/home/toki/miniforge3/envs/hifiasm/bin:/home/toki/miniforge3/envs/jf/bin:/home/toki/miniforge3/envs/sinederella/bin:/home/toki/miniforge3/envs/asm25/bin:/home/toki/miniforge3/envs/lrtools/bin:$PATH
export TMPDIR=/home/toki/tmp
T=${T:-32}
if [ "$MODE" = toy ]; then ASMS="dep1"; else ASMS="dep1 dep2"; fi
# NEW_ASM = none: only the hifiasm 0.16.1 -l0 sources (first in source order: wgsp)
if [ "$NEW" = none ]; then SRC="wgsp=$WP wgsa=$WA"; FIRST=wgsp; else SRC="hicanu=$NEW wgsp=$WP wgsa=$WA"; FIRST=hicanu; fi

for a in $ASMS; do
  if [ $a = dep1 ]; then GF=$D1; CA=/home/toki/gof/gf_all/flankcn/gap_causes.tsv; else GF=$D2; CA=$R1/causes2/gap_causes.tsv; fi
  echo "#### gapfill $a $(date '+%F %T')"
  python3 $W/gapfill.py $GF gf_$a $SRC --threads $T
  echo "#### fill_cn $a $(date '+%F %T')"
  meryl-lookup -wig-count -sequence gf_$a/fills.fa -mers /home/toki/gof/sqc/hifi.meryl -output gf_$a/fills.wig > gf_$a/meryl.log 2>&1
  python3 $C/fill_cn.py gf_$a/fills.fa gf_$a/fills.wig /home/toki/gof/ends/kmer/gc_calib.txt > gf_$a/fill_cn.tsv
  rm -f gf_$a/fills.wig
  echo "#### measure $a $(date '+%F %T')"
  python3 $C/asm2_measure.py gf_$a/gaps.tsv gf_$a/fill_cn.tsv gf_$a/fills.fa $CA $R1/repairs.tsv $a $a
  G2=$(cat $a.gaps2.txt)
  if [ "$MODE" = toy ]; then G2=${TOY_GAPS:-$(cut -d, -f1-6 $a.gaps2.txt)}; fi
  [ -n "$G2" ] || { echo "no gaps for round 2 in $a"; continue; }
  echo "#### plan $a $(date '+%F %T')"
  python3 $C/repair_gaps.py plan $GF gf_$a/gaps.tsv gf_$a/fills.fa gf_$a/fills.fa gf_$a/fill_cn.tsv $a.causes2.tsv $FIRST r2_$a --gaps "$G2"
  n=$(grep -c '>' r2_$a/cand.fa || true)
  [ "$n" -gt 0 ] || { echo "nothing built in $a"; continue; }
  jellyfish count -m 25 -C -s $((4 * n + 1000)) -t 16 --if r2_$a/cand.fa -o r2_$a/cand.jf $GF r2_$a/pieces.fa
  jellyfish dump -c r2_$a/cand.jf > r2_$a/cand.counts
  python3 $C/repair_gaps.py words r2_$a r2_$a/cand.counts
  rm -f r2_$a/cand.jf r2_$a/cand.counts
done

echo "#### fish $(date '+%F %T')"
if [ "$MODE" = toy ]; then
  pigz -dc $R1/fished.fq.gz > fished.fq
else
  cat $(for a in $ASMS; do [ -s r2_$a/anchors.txt ] && echo r2_$a/anchors.txt; done) | sort -u > anchors.txt
  wc -l anchors.txt
  pigz -dc -p 8 $R | paste - - - - | /home/toki/gof/scfish/scfish anchors.txt 3 \
    | awk -F'\t' '{split($1,h," "); print "@" substr(h[1],2) "\n" $2 "\n+\n" $4}' > fished.fq
fi
echo "fished reads: $(( $(wc -l < fished.fq) / 4 ))"

for a in $ASMS; do
  [ -s r2_$a/junctions.tsv ] || continue
  if [ $a = dep1 ]; then GF=$D1; else GF=$D2; fi
  echo "#### verify $a $(date '+%F %T')"
  python3 $C/repair_gaps.py verify $GF r2_$a fished.fq $a $T --edited-only
  rm -f r2_$a/$a.repaired.fa r2_$a/$a.genesearch.fa       # round-2-only sequences; merged below
  echo "#### merge $a $(date '+%F %T')"
  python3 $C/asm2_merge.py $GF $R1/$a r2_$a $a . | tee $a.merge.txt
done
pigz -p 8 -f fished.fq
echo "#### done $(date '+%F %T')"
