#!/bin/bash
set -euo pipefail
C=/home/toki/gof/code
W=/home/toki/gof/toy
mkdir -p $W
cd $W || exit 1
export PATH=/home/toki/toga2/env_mlc2/bin:/home/toki/miniforge3/envs/sinederella/bin:$PATH
export TMPDIR=/home/toki/tmp
T=/home/toki/toga2/b2t
TO=$T/full/toga_out
python3 $C/tests/make_toy.py $TO $T/sbe_ref_isoforms.tsv /home/toki/toga2/inputs/trizona_nodot.2bit $W 100
# betulina contigs that carry the toy genes, for the toy whole-genome alignment (edge scale)
python3 - "$T/sbe_ref_isoforms.tsv" "$T/sbe_ref_annotation.bed" genes.txt > toy_contigs.txt <<'PY'
import sys
iso, bed, genes = sys.argv[1:4]
want = set(l.split()[0] for l in open(genes))
tx = {l.split("\t")[1].strip(): l.split("\t")[0] for l in open(iso)}
print("\n".join(sorted({l.split("\t")[0] for l in open(bed) if tx.get(l.split("\t")[3]) in want})))
PY
twoBitToFa -seqList=toy_contigs.txt $T/sbe.2bit $TMPDIR/toy_sbe.fa
for v in keep del gap; do
  faToTwoBit $v.fa $v.2bit
  minimap2 -x asm20 -d $v.mmi $v.fa 2> $v.idx.log
  minimap2 -c -x asm20 -t 16 --secondary=no $v.mmi $TMPDIR/toy_sbe.fa > wga_$v.paf 2> wga_$v.log
  python3 $C/gof.py --toga $TO --ref-bed $T/sbe_ref_annotation.bed --isoforms $T/sbe_ref_isoforms.tsv \
    --ref-2bit $T/sbe.2bit --query-2bit $v.2bit --query-index $v.mmi --genes genes.txt --out $W/out_$v \
    --threads 8 --no-neighbours --wga-paf wga_$v.paf
  echo "== $v: f1000 / f10000 / edge verdicts"
  awk -F'\t' 'NR>1{print $5"\t"$8"\t"($12==""?"-":$12)}' out_$v/verdicts.tsv | sort | uniq -c
done
rm -f $TMPDIR/toy_sbe.fa
for v in keep del gap; do echo "== $v: gene call"; awk -F'\t' 'NR>1{print $NF}' out_$v/verdicts.tsv | sort | uniq -c; done
