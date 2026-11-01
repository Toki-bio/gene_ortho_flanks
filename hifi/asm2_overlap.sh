#!/bin/bash
# Round 2 addition: overlap merges at N runs whose contig ends carry the same sequence (asm2_overlap.py), verified
# by reads crossing the merged junction (repair_gaps.py words/verify, unchanged rule), merged with round 1 and the
# round-2 new-source repairs (asm2_merge.py).
#  dep1: the 17 'tandem' gaps with read-measured content <= -2 kb (recheck_label.py relabel.tsv), searched around E.
#  dep1 blind: the same 17 without E (validation: does the sequence alone find the read-measured overlap?).
#  dep1 fallback: no exact duplicate -> array over-length trim of E bp (asm2_overlap.py docstring).
#  dep2: all 'tandem' gaps (causes2), blind exact-duplicate search only (no read-measured content for dep2),
#  kept only at overlap identity >= 0.8.
# Usage: asm2_overlap.sh OUT_DIR [toy]   (toy: dep1, 3 gaps, reads = recheck fished reads)
set -euo pipefail
O=$1
MODE=${2:-full}
mkdir -p "$O"
cd "$O" || exit 1
C=/home/toki/gof/code/hifi
R=/home/toki/gof/reads/ERR16751478.fastq.gz
R1=/home/toki/gof/repair
R2=/home/toki/gof/asm2/r2_wgs
REL=/home/toki/gof/recheck/main/res/relabel.tsv
D1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
D2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
export PATH=/home/toki/miniforge3/envs/hifiasm/bin:/home/toki/miniforge3/envs/jf/bin:/home/toki/miniforge3/envs/lrtools/bin:$PATH
export TMPDIR=/home/toki/tmp
T=${T:-16}
S1=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["old_label"]=="tandem" && $h["content_median"]!="" && $h["content_median"]+0 <= -2000 {printf "%s%s:%d:%d:%d", (n++?",":""), $1, -$h["content_median"], $h["array_L"], $h["array_R"]}' $REL)
if [ "$MODE" = toy ]; then S1=$(echo $S1 | cut -d, -f4,6,8); fi
echo "dep1 gaps: $S1"
python3 $C/asm2_overlap.py $D1 /home/toki/gof/gf_all/dep1/gaps.tsv ov_dep1 "$S1"
python3 $C/asm2_overlap.py $D1 /home/toki/gof/gf_all/dep1/gaps.tsv ov_dep1_blind "$(echo $S1 | sed 's/:[0-9:]*//g')" > /dev/null
ASMS=dep1
if [ "$MODE" != toy ]; then
  S2=$(awk -F'\t' 'NR>1 && $NF=="tandem"{printf "%s%s", (n++?",":""), $1}' $R1/causes2/gap_causes.tsv)
  echo "dep2 tandem gaps: $(echo $S2 | tr ',' '\n' | wc -l)"
  python3 $C/asm2_overlap.py $D2 /home/toki/gof/gf_all/dep2/gaps.tsv ov_dep2 "$S2" > /dev/null
  ASMS="dep1 dep2"
fi
for a in $ASMS; do
  if [ $a = dep1 ]; then GF=$D1; else GF=$D2; fi
  n=$(grep -c '>' ov_$a/cand.fa || true)
  [ "$n" -gt 0 ] || continue
  jellyfish count -m 25 -C -s $((4 * n + 1000)) -t 16 --if ov_$a/cand.fa -o ov_$a/cand.jf $GF
  jellyfish dump -c ov_$a/cand.jf > ov_$a/cand.counts
  python3 $C/repair_gaps.py words ov_$a ov_$a/cand.counts
  rm -f ov_$a/cand.jf ov_$a/cand.counts
done
echo "#### fish $(date '+%F %T')"
if [ "$MODE" = toy ]; then
  ln -sf /home/toki/gof/recheck/main/fished.fq fished.fq
else
  cat $(for a in $ASMS; do [ -s ov_$a/anchors.txt ] && echo ov_$a/anchors.txt; done) | sort -u > anchors.txt
  pigz -dc -p 8 $R | paste - - - - | /home/toki/gof/scfish/scfish anchors.txt 3 \
    | awk -F'\t' '{split($1,h," "); print "@" substr(h[1],2) "\n" $2 "\n+\n" $4}' > fished.fq
fi
echo "fished reads: $(( $(wc -l < fished.fq) / 4 ))"
for a in $ASMS; do
  [ -s ov_$a/junctions.tsv ] || continue
  if [ $a = dep1 ]; then GF=$D1; else GF=$D2; fi
  echo "#### verify $a $(date '+%F %T')"
  python3 $C/repair_gaps.py verify $GF ov_$a fished.fq $a $T --edited-only
  rm -f ov_$a/$a.repaired.fa ov_$a/$a.genesearch.fa
  echo "#### merge $a $(date '+%F %T')"
  python3 $C/asm2_merge.py $GF $R1/$a $R2/r2_$a,ov_$a $a . | tee $a.merge.txt
done
[ -L fished.fq ] || pigz -p 8 -f fished.fq
echo "#### done $(date '+%F %T')"
