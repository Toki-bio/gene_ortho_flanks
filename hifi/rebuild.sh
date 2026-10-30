#!/bin/bash
# Rebuild query regions (gaps, deletion sites, controls) from the query's own HiFi reads,
# sineclose-style: unique 25-mer baits -> scfish -> per-region hifiasm -> evaluate.
# Usage: rebuild.sh STAGE   (baits | fish | assemble | evaluate | all)
set -euo pipefail
W=${W:-/home/toki/gof/hifi}
cd $W || exit 1
[ "$(pwd)" = "$W" ] || exit 1
export PATH=/home/toki/toga2/env_mlc2/bin:/home/toki/miniforge3/envs/sinederella/bin:/home/toki/miniforge3/envs/hifiasm/bin:$PATH
export TMPDIR=/home/toki/tmp
Q2=/home/toki/toga2/inputs/trizona_nodot.2bit
QFA=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
R2=/home/toki/toga2/b2t/sbe.2bit
RFA=/home/toki/Sicista2026/assembly/primary_eval/hifiasm_primary.p_ctg.fa
READS=${READS:-/home/toki/gof/reads/ERR16751478.fastq.gz}
SCFISH=/home/toki/gof/scfish/scfish
PAD=10000
K=25
MIN_HITS=3
step() { echo "#### $* -- $(date '+%F %T')"; }
stage=${1:-all}

if [ $stage = baits ] || [ $stage = all ]; then
  step baits
  # query side: region +- PAD; reference side: the reference interval
  awk -F'\t' -v P=$PAD 'NR>1{s=$4-P; if(s<0)s=0; print $3"\t"s"\t"$5+P"\t"$1"|q"; print $6"\t"$7"\t"$8"\t"$1"|r"}' regions.tsv > bait.bed
  awk -F'\t' '$4 ~ /\|q$/' bait.bed > bait_q.bed
  awk -F'\t' '$4 ~ /\|r$/' bait.bed > bait_r.bed
  twoBitToFa -bed=bait_q.bed $Q2 bait_q.fa
  twoBitToFa -bed=bait_r.bed $R2 bait_r.fa
  cat bait_q.fa bait_r.fa > bait.fa
  # genome counts of the bait k-mers only; keep k-mers that occur at most once in each genome
  jellyfish count -m $K -C -s 100M -t 16 --if bait.fa -o q.jf $QFA
  jellyfish count -m $K -C -s 100M -t 16 --if bait.fa -o r.jf $RFA
  jellyfish dump -c q.jf > q_counts.txt
  jellyfish dump -c r.jf > r_counts.txt
  python3 - <<'PY'
import re
K = 25
comp = str.maketrans("ACGT", "TGCA")
def canon(k):
    r = k.translate(comp)[::-1]
    return min(k, r)
q = dict(l.split() for l in open("q_counts.txt"))
r = dict(l.split() for l in open("r_counts.txt"))
seqs, name = {}, None
for l in open("bait.fa"):
    l = l.strip()
    if l.startswith(">"):
        name = l[1:]; seqs[name] = []
    else:
        seqs[name].append(l.upper())
keep, total = set(), 0
for n, parts in seqs.items():
    s = "".join(parts)
    for i in range(len(s) - K + 1):
        k = s[i:i + K]
        if not re.fullmatch("[ACGT]+", k):
            continue
        total += 1
        c = canon(k)
        if int(q.get(c, 0)) <= 1 and int(r.get(c, 0)) <= 1:
            keep.add(c)
open("anchors.txt", "w").write("\n".join(sorted(keep)) + "\n")
print(f"bait k-mer positions {total}, unique anchors kept {len(keep)}")
PY
  rm -f q.jf r.jf
fi

if [ $stage = fish ] || [ $stage = all ]; then
  step fish
  [ -s $READS ] || { echo "reads not downloaded"; exit 1; }
  UNZ="gzip -dc"; command -v pigz > /dev/null && UNZ="pigz -dc -p 4"
  $UNZ $READS | paste - - - - | $SCFISH anchors.txt $MIN_HITS > fished.tsv
  awk -F'\t' '{print $1"\n"$2"\n+\n"$4}' fished.tsv > fished.fq
  echo "fished reads $(wc -l < fished.tsv), bases $(awk -F'\t' '{s+=length($2)} END{print s}' fished.tsv)"
fi

if [ $stage = assemble ] || [ $stage = all ]; then
  step assemble
  # each read goes to the region(s) whose bait it maps to
  minimap2 -x map-hifi -t 16 bait.fa fished.fq 2> map.log | awk -F'\t' '$12>=0 && ($4-$3)>=1000 {split($6,a,"|"); print $1"\t"a[1]}' | sort -u > read2region.tsv
  rm -rf asm && mkdir asm
  for id in $(cut -f2 read2region.tsv | sort -u); do
    mkdir -p asm/$id
    awk -v id=$id '$2==id{print $1}' read2region.tsv | sort -u > asm/$id/names.txt
    seqkit grep -f asm/$id/names.txt fished.fq > asm/$id/reads.fq 2>/dev/null
    ( cd asm/$id || exit 1
      hifiasm -o h -t 8 -l0 reads.fq > hifiasm.log 2>&1 || true
      [ -s h.bp.p_ctg.gfa ] && awk '$1=="S"{print ">"$2"\n"$3}' h.bp.p_ctg.gfa > contigs.fa || : > contigs.fa
      echo "$id reads $(wc -l < names.txt) contigs $(grep -c '>' contigs.fa)" )
  done
fi

if [ $stage = evaluate ] || [ $stage = all ]; then
  step evaluate
  python3 /home/toki/gof/code/hifi/evaluate.py regions.tsv bait.fa asm $PAD > rebuild_results.tsv
  column -t rebuild_results.tsv
fi
step finished
