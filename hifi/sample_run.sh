#!/bin/bash
# One-copy vs two-copy fills: sample gaps (sample_ends.py) -> scfish pass over the HiFi reads -> classify_ends.py
set -euo pipefail
O=/home/toki/gof/gf_all/sample
mkdir -p $O
cd $O || exit 1
C=/home/toki/gof/code/hifi
G=/home/toki/gof/gf_all
export PATH=/home/toki/miniforge3/envs/hifiasm/bin:/home/toki/miniforge3/envs/jf/bin:/home/toki/miniforge3/envs/sinederella/bin:/home/toki/miniforge3/envs/lrtools/bin:$PATH
echo "#### sample $(date '+%F %T')"
python3 $C/sample_ends.py /home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna $G/dep1/gaps.tsv $G/fill_cn.tsv dep2 $G/dep1/fills.fa $O 25 15
echo "#### fish $(date '+%F %T')"
pigz -dc -p 8 /home/toki/gof/reads/ERR16751478.fastq.gz | paste - - - - | /home/toki/gof/scfish/scfish anchors.txt 3 \
  | awk -F'\t' '{split($1,h," "); print "@" substr(h[1],2) "\n" $2 "\n+\n" $4}' > fished.fq
echo "#### classify $(date '+%F %T')"
python3 $C/classify_ends.py fills.fa fished.fq cls 32
echo "#### done $(date '+%F %T')"
