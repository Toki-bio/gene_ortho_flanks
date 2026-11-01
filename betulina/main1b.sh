#!/bin/bash
# Rest of main1.sh (main1 broke at its trf line because the script file was overwritten while it ran).
set -euo pipefail
cd /home/toki/gof/betulina/gw || exit 1
E=/home/toki/miniforge3/envs
export PATH=$E/sinederella/bin:$PATH
MM=$E/hifiasm/bin/minimap2
C=/home/toki/gof/betulina/code
HC=/home/toki/gof/code/hifi
P=/home/toki/Sicista2026/assembly/primary_eval/hifiasm_primary.p_ctg.fa
H1=/home/toki/Sicista2026/assembly/hifiasm_out/sicista_noscaf.p_ctg.fa
T1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
TGAPS=/home/toki/gof/gf_all/dep1/gaps.tsv
TFL=/home/toki/gof/gf_all/flankcn/flanks.fa
date
# TRF on the end flanks, 8 parts in parallel
mkdir -p trfp
awk -v n=8 '/^>/{i++} {print > ("trfp/p" (i % n) ".fa")}' ends.ends.fa
for f in trfp/p*.fa; do ( cd trfp && trf $(basename $f) 2 7 7 80 10 50 500 -h -ngs > $(basename $f .fa).dat || true ) & done
wait
cat trfp/p*.dat > ends.trf.dat
echo "trf seqs with hits: $(grep -c '^@' ends.trf.dat)"
python3 $C/random_sites.py $T1.fai $TGAPS 1000 7 > rnd.tsv
python3 $HC/gap_flank_cn.py prep $T1 rnd.tsv rnd_flanks.fa
cat $TFL rnd_flanks.fa ends.ends.fa > q_all.fa
$MM -x asm20 -t 12 $P q_all.fa 2> qa_p.mm.log > qa_vs_p.paf
$MM -x asm20 -t 12 $H1 q_all.fa 2> qa_h1.mm.log > qa_vs_h1.paf
$MM -x asm20 -t 12 bet_hap2.fa q_all.fa 2> qa_h2.mm.log > qa_vs_h2.paf
$MM -x asm20 -t 12 $T1 ends.ends.fa 2> ends_t1.mm.log > ends_vs_t1.paf
date; echo F_done
