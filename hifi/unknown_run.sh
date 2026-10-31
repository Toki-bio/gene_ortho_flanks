#!/bin/bash
# Characterise unspanned gap flanks vs one/two-copy fills vs random controls: RepeatMasker (rodent Dfam), TRF, k-mer CN.
set -euo pipefail
O=/home/toki/gof/gf_all/unk
mkdir -p $O
cd $O || exit 1
C=/home/toki/gof/code/hifi
G=/home/toki/gof/gf_all
export PATH=/home/toki/miniforge3/envs/sinederella/bin:/home/toki/miniforge3/envs/asm25/bin:$PATH
export TMPDIR=/home/toki/tmp
D1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
echo "#### prep $(date '+%F %T')"
python3 $C/unknown_prep.py $D1 $G/dep1/gaps.tsv $G/dep1/fills.fa $G/fill_cn.tsv dep2 /home/toki/refs/dfam_current/curated_consensi.fa /home/toki/refs/dfam_current/curated_meta.tsv $O
grep -c '>' feat.fa
echo "#### trf $(date '+%F %T')"
[ -s feat.trf.dat ] || trf feat.fa 2 7 7 80 10 50 2000 -h -d -ngs > feat.trf.dat || true
echo "#### repeatmasker $(date '+%F %T')"
RepeatMasker -lib rodent_lib.fa -pa 16 -dir rm -no_is feat.fa > rm.log 2>&1 || echo "RepeatMasker failed (see rm.log); repeat classes left empty"
[ -s rm/feat.fa.out ] || { mkdir -p rm; : > rm/feat.fa.out; }
echo "#### meryl $(date '+%F %T')"
meryl-lookup -wig-count -sequence feat.fa -mers /home/toki/gof/sqc/hifi.meryl -output feat.wig > meryl.log 2>&1
echo "#### summary $(date '+%F %T')"
python3 $C/unknown_feat.py feat.fa rm/feat.fa.out feat.trf.dat feat.wig $O > summary.txt
cat summary.txt
echo "#### done $(date '+%F %T')"
