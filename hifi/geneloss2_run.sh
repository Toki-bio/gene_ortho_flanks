#!/bin/bash
# geneloss_repaired.py full run: 191 genes (b2t deletion/conflict/assembly_gap + 20 FI controls) on deposited hap1/hap2
# before repair (orig1/orig2), repaired (dep1/dep2) and gene-search (dep1gs/dep2gs) versions.
set -euo pipefail
cd /home/toki/gof/geneloss2 || exit 1
export PATH=/home/toki/miniforge3/envs/hifiasm/bin:/home/toki/toga2/env_mlc2/bin:/home/toki/miniforge3/envs/sinederella/bin:$PATH
export TMPDIR=/home/toki/tmp
C=/home/toki/gof/code
G=/home/toki/gof/gf_all
R=/home/toki/gof/repair
D1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
D2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
O=prep
for a in orig1 orig2 dep1 dep2 dep1gs dep2gs; do
  case $a in
    orig1) F=$D1; OR=$D1; GP=$G/dep1/gaps.tsv; CA=$G/flankcn/gap_causes.tsv;;
    orig2) F=$D2; OR=$D2; GP=$G/dep2/gaps.tsv; CA=$R/causes2/gap_causes.tsv;;
    dep1) F=$R/dep1.repaired.fa; OR=$D1; GP=$G/dep1/gaps.tsv; CA=$G/flankcn/gap_causes.tsv;;
    dep2) F=$R/dep2.repaired.fa; OR=$D2; GP=$G/dep2/gaps.tsv; CA=$R/causes2/gap_causes.tsv;;
    dep1gs) F=$R/dep1.genesearch.fa; OR=$D1; GP=$G/dep1/gaps.tsv; CA=$G/flankcn/gap_causes.tsv;;
    dep2gs) F=$R/dep2.genesearch.fa; OR=$D2; GP=$G/dep2/gaps.tsv; CA=$R/causes2/gap_causes.tsv;;
  esac
  echo "#### $a $(date '+%F %T')"
  minimap2 -c --cs -x asm20 -N 5 -t 32 $F $O/query.fa > $O/$a.paf 2> $O/$a.mm.log
  python3 $C/hifi/geneloss_repaired.py call $O $a $O/$a.paf $F $OR $GP $CA
done
python3 $C/hifi/geneloss_repaired.py combine $O
cp $O/final.tsv $O/summary.txt .
echo "#### done $(date '+%F %T')"
