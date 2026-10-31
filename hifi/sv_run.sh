#!/bin/bash
# Step 1b: windows + words (sv_prep.py) -> one scfish pass over the HiFi reads -> sv_window.py
set -euo pipefail
O=/home/toki/gof/ends/sv
mkdir -p $O
cd $O || exit 1
C=/home/toki/gof/code/hifi
P=/home/toki/gof/hifi/paralog
A=/home/toki/gof/asm25
export PATH=/home/toki/miniforge3/envs/hifiasm/bin:/home/toki/miniforge3/envs/jf/bin:/home/toki/miniforge3/envs/sinederella/bin:/home/toki/miniforge3/envs/lrtools/bin:$PATH
for t in minimap2 jellyfish samtools pigz; do command -v $t > /dev/null || { echo "missing $t"; exit 1; }; done
echo "#### prep $(date '+%F %T')"
python3 $C/sv_prep.py $O hap1=$A/tri.hic.hap1.p_ctg.fa:$P/fills_vs_asm25hap1.paf hap2=$A/tri.hic.hap2.p_ctg.fa:$P/fills_vs_asm25hap2.paf
echo "#### fish $(date '+%F %T')"
pigz -dc -p 8 /home/toki/gof/reads/ERR16751478.fastq.gz | paste - - - - | /home/toki/gof/scfish/scfish anchors.txt 3 \
  | awk -F'\t' '{split($1,h," "); print ">" substr(h[1],2) "\n" $2}' > fished.fa
grep -c '>' fished.fa
echo "#### windows $(date '+%F %T')"
python3 $C/sv_window.py windows.fa fished.fa out 32 > sv_window.txt
echo "#### done $(date '+%F %T')"
