#!/bin/bash
# After recheck_reads.sh: re-label from read evidence (recheck_label.py) and locate the other copies of high-copy
# flanks + where read tails leaving them go (recheck_copies.py).
# Usage: recheck_post.sh DIR
set -euo pipefail
D=$1
cd "$D" || exit 1
C=/home/toki/gof/code/hifi
export PATH=/home/toki/miniforge3/envs/sinederella/bin:$PATH
[ -s res/relabel.tsv ] || python3 $C/recheck_label.py win.fa h1.paf h1.trf.dat res targets.tsv > res/relabel.tsv
python3 $C/recheck_copies.py prep win.fa /home/toki/gof/gf_all/flankcn/gap_causes.tsv high_copy,tandem copies_segs.fa 3000
[ -s copies.paf ] || minimap2 -c -x asm20 -N 200 -p 0.1 -t 16 /home/toki/gof/trizona.asm20.mmi copies_segs.fa res/tails.fa > copies.paf 2> copies_mm.log
grep '>' copies_segs.fa | cut -c2- > copies_segs.ids
awk -F'\t' 'NR==FNR{s[$1]=1;next} ($1 in s)' copies_segs.ids copies.paf > copies_segs.paf
awk -F'\t' 'NR==FNR{s[$1]=1;next} !($1 in s)' copies_segs.ids copies.paf > copies_tails.paf
python3 $C/recheck_copies.py table copies_segs.paf copies_tails.paf /home/toki/gof/gf_all/flankcn/gap_causes.tsv high_copy,tandem res/tails.fa > res/copies.tsv
cut -f2,20 res/relabel.tsv | sort | uniq -c
