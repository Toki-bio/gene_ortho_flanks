#!/bin/bash
# Two-copy proof for repaired gap pieces (prove_two_copy.py):
#   prep -> HiFi anchors (unique 25-mers of +-200 kb regions) -> fish HiFi -> genome-wide map to both gene-search
#   assemblies -> het sites from reads -> test A (phased reads crossing the junctions)
#   -> Hi-C allele 31-mers -> fish Hi-C pairs -> map mates single-end -> test B (allele-specific contacts)
#   -> test C (betulina primary carries piece +-10 kb collinearly)
# Usage: prove_run.sh OUT_DIR [toy]
set -euo pipefail
O=$1
MODE=${2:-full}
mkdir -p "$O"
cd "$O" || exit 1
C=/home/toki/gof/code/hifi
P="python3 $C/prove_two_copy.py"
R=/home/toki/gof/repair
G1=$R/dep1.genesearch.fa
G2=$R/dep2.genesearch.fa
HIFI=/home/toki/gof/reads/ERR16751478.fastq.gz
HC1=/home/toki/gof/hic/ERR18010482_1.fastq.gz
HC2=/home/toki/gof/hic/ERR18010482_2.fastq.gz
BET=/home/toki/Sicista2026/assembly/primary_eval/hifiasm_primary.p_ctg.fa
IDX=/home/toki/gof/proof2copy/idx
E=/home/toki/miniforge3/envs
export PATH=$E/hifiasm/bin:$E/jf/bin:$E/lrtools/bin:$PATH
T=32
HETCTL="--hetctl TMEM114:/home/toki/gof/views_het/tmem_t2.fa:15708:79274 --hetctl SMCO1:/home/toki/gof/views_het/smco_t2.fa:15502:23296"
mkdir -p $IDX

echo "#### prep $(date '+%F %T')"
if [ "$MODE" = toy ]; then
  mkdir -p pre
  $P prep $G1 $R/dep1.genesearch.bed $G2 $R/dep2.genesearch.bed pre $HETCTL
  ONLY=$(awk -F'\t' 'NR>1 && !seen[$2]++ {printf "%s,", $1}' pre/loci.tsv)
  $P prep $G1 $R/dep1.genesearch.bed $G2 $R/dep2.genesearch.bed . $HETCTL --only "${ONLY%,}"
  HIFI_CMD="zcat $R/fished.fq.gz"   # reads already fished around every repair junction
  HIC_N=16000000
else
  $P prep $G1 $R/dep1.genesearch.bed $G2 $R/dep2.genesearch.bed . $HETCTL
  HIFI_CMD="pigz -dc -p 8 $HIFI"
  HIC_N=0
fi
cut -f2 loci.tsv | sort | uniq -c

echo "#### indexes $(date '+%F %T')"
for a in 1 2; do
  G=$([ $a = 1 ] && echo $G1 || echo $G2)
  [ -s $IDX/dep${a}g.hifi.mmi ] || minimap2 -x map-hifi -t 8 -d $IDX/dep${a}g.hifi.mmi $G 2> $IDX/dep${a}g.hifi.log
  [ -s $IDX/dep${a}g.sr.mmi ] || minimap2 -x sr -t 8 -d $IDX/dep${a}g.sr.mmi $G 2> $IDX/dep${a}g.sr.log
done

echo "#### HiFi anchors $(date '+%F %T')"
n=$(grep -c '>' hifi_cand.fa)
jellyfish count -m 25 -C -s $((4 * n + 1000)) -t 16 --if hifi_cand.fa -o hifi.jf $G1 $G2
jellyfish dump -c hifi.jf > hifi.counts
$P anchors hifi_cand.fa hifi.counts 1 2 > hifi_anchors.txt
rm -f hifi.jf hifi.counts
wc -l < hifi_anchors.txt

echo "#### fish HiFi $(date '+%F %T')"
eval "$HIFI_CMD" | paste - - - - | /home/toki/gof/scfish/scfish hifi_anchors.txt 3 \
  | awk -F'\t' '{split($1,h," "); print "@" substr(h[1],2) "\n" $2 "\n+\n" $4}' | gzip -1 > hifi_fished.fq.gz
echo "fished HiFi reads: $(( $(zcat hifi_fished.fq.gz | wc -l) / 4 ))"

echo "#### map HiFi $(date '+%F %T')"
for a in 1 2; do
  minimap2 -x map-hifi -c --cs --secondary=no -t $T $IDX/dep${a}g.hifi.mmi hifi_fished.fq.gz 2> map_hifi_dep$a.log \
    | $P filter regions.bed dep$a | gzip -1 > hifi_dep$a.paf.gz
  echo "dep$a alignments kept: $(zcat hifi_dep$a.paf.gz | wc -l)"
done

echo "#### test A $(date '+%F %T')"
$P het loci.tsv regions.bed hifi_dep1.paf.gz hifi_dep2.paf.gz > het.tsv
echo "het sites: $(( $(wc -l < het.tsv) - 1 ))"
$P phase loci.tsv het.tsv hifi_dep1.paf.gz hifi_dep2.paf.gz testA
cat testA.summary.txt

echo "#### Hi-C anchors $(date '+%F %T')"
$P hicanc loci.tsv het.tsv regions.fa .
n=$(grep -c '>' hic_cand.fa)
jellyfish count -m 31 -C -s $((4 * n + 1000)) -t 16 --if hic_cand.fa -o hic.jf $G1 $G2
jellyfish dump -c hic.jf > hic.counts
$P anchors hic_cand.fa hic.counts 0 2 > hic_anchors.txt
rm -f hic.jf hic.counts
echo "Hi-C sites: $(wc -l < hic_sites.tsv), anchors: $(wc -l < hic_anchors.txt)"

echo "#### fish Hi-C $(date '+%F %T')"
if [ "$HIC_N" -gt 0 ]; then
  paste <((pigz -dc -p 4 $HC1 || true) | head -n $HIC_N | paste - - - -) <((pigz -dc -p 4 $HC2 || true) | head -n $HIC_N | paste - - - -) \
    | /home/toki/gof/scfish/scfish hic_anchors.txt 1 > hic_fished.tsv
else
  paste <(pigz -dc -p 4 $HC1 | paste - - - -) <(pigz -dc -p 4 $HC2 | paste - - - -) \
    | /home/toki/gof/scfish/scfish hic_anchors.txt 1 > hic_fished.tsv
fi
awk -F'\t' '{split($1,h," "); n=substr(h[1],2); print "@" n "_1\n" $2 "\n+\n" $4 "\n@" n "_2\n" $6 "\n+\n" $8}' hic_fished.tsv \
  | gzip -1 > hic_fished.fq.gz
echo "fished Hi-C pairs: $(wc -l < hic_fished.tsv)"
rm -f hic_fished.tsv

echo "#### map Hi-C $(date '+%F %T')"
for a in 1 2; do
  minimap2 -x sr -c --cs --secondary=no -t $T $IDX/dep${a}g.sr.mmi hic_fished.fq.gz 2> map_hic_dep$a.log | gzip -1 > hic_dep$a.paf.gz
done
$P hic loci.tsv hic_sites.tsv hic_dep1.paf.gz hic_dep2.paf.gz testB
cat testB.summary.txt

echo "#### test C $(date '+%F %T')"
$P betq loci.tsv $G1 $G2 500000 > bet_q.fa
minimap2 -x asm20 -c --secondary=no -t $T $BET bet_q.fa 2> map_bet.log > bet.paf
$P bet loci.tsv bet.paf > testC.tsv
awk -F'\t' 'NR>1{print $2, $NF}' testC.tsv | sort | uniq -c
echo "#### done $(date '+%F %T')"
