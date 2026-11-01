#!/bin/bash
# Re-examination of weak gap labels: windows around the selected gaps (+ alt fill windows) -> unique 25-mers over
# hap1+hap2 -> one scfish pass over the HiFi reads (toy: first 400k reads).
# Usage: recheck_fish.sh OUT_DIR [toy]
set -euo pipefail
O=$1
MODE=${2:-full}
mkdir -p "$O"
cd "$O" || exit 1
C=/home/toki/gof/code/hifi
G=/home/toki/gof/gf_all
R=/home/toki/gof/reads/ERR16751478.fastq.gz
D1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
export DEP2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
export ASM_H1=/home/toki/gof/asm25/tri.hic.hap1.p_ctg.fa
export ASM_H2=/home/toki/gof/asm25/tri.hic.hap2.p_ctg.fa
export PATH=/home/toki/miniforge3/envs/hifiasm/bin:/home/toki/miniforge3/envs/jf/bin:/home/toki/miniforge3/envs/sinederella/bin:/home/toki/miniforge3/envs/lrtools/bin:$PATH
if [ "$MODE" = toy ]; then
  SEL=none; EXTRA=gap0002,gap0011
else
  SEL=join,tandem,high_copy,ga_tract,fill_repeat
  # one-copy examples: one_copy_fill gaps whose dep2 fill is 3-12 kb
  EXTRA=$(awk -F'\t' 'NR==FNR{if($14=="one_copy_fill")ok[$1]=1;next} ($2=="dep2" && ($1 in ok) && $3>=3000 && $3<=12000){print $1}' \
    $G/flankcn/gap_causes.tsv $G/fill_cn.tsv | head -14 | paste -sd,)
fi
echo "#### prep $(date '+%F %T') extra=$EXTRA"
python3 $C/recheck_prep.py prep $D1 $G/dep1/gaps.tsv $G/flankcn/gap_causes.tsv $SEL $EXTRA . 40000 15
n=$(grep -c '>' cand.fa)
echo "#### count $(date '+%F %T') $n words"
jellyfish count -m 25 -C -s $((2 * n + 1000)) -t 16 --if cand.fa -o cand.jf $D1 $DEP2
jellyfish dump -c cand.jf > cand.counts
python3 $C/recheck_prep.py words . cand.counts 2
rm -f cand.jf cand.counts
echo "#### fish $(date '+%F %T')"
if [ "$MODE" = toy ]; then
  pigz -dc -p 4 $R | head -1600000 | paste - - - - | /home/toki/gof/scfish/scfish anchors.txt 3 \
    | awk -F'\t' '{split($1,h," "); print "@" substr(h[1],2) "\n" $2 "\n+\n" $4}' > fished.fq || true
else
  pigz -dc -p 8 $R | paste - - - - | /home/toki/gof/scfish/scfish anchors.txt 3 \
    | awk -F'\t' '{split($1,h," "); print "@" substr(h[1],2) "\n" $2 "\n+\n" $4}' > fished.fq
fi
echo "fished reads: $(( $(wc -l < fished.fq) / 4 ))"
echo "#### done $(date '+%F %T')"
