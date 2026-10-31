#!/bin/bash
# Hi-C scaffolding (YaHS) and QC of one haplotype assembly.
# chromap Hi-C mapping streamed into a name-sorted BAM (no SAM on disk), YaHS, gfastats,
# compleasm, Merqury QV/completeness against a meryl db of the HiFi reads.
# COMPLEASM_LIB: offline library dir (busco-data.ezlab.org is unreachable from therioserver) holding
# file_versions.tsv, mammalia_odb10/ (BUSCO lineage) and the .done markers compleasm checks.
# Usage: scaffold_qc.sh CONTIGS.fa HIC_R1 HIC_R2 OUT_DIR COMPLEASM_LIB MERYL_DB [THREADS]
set -euo pipefail
A=$(readlink -f $1); R1=$(readlink -f $2); R2=$(readlink -f $3); O=$4
LIN=$(readlink -f $5); MDB=$(readlink -f $6); T=${7:-48}
mkdir -p $O
cd $O || exit 1
export PATH=/home/toki/miniforge3/envs/asm25/bin:$PATH
export TMPDIR=/home/toki/tmp
export MERQURY=$(ls -d /home/toki/miniforge3/envs/asm25/share/merqury*)
step() { echo "#### $* -- $(date '+%F %T')"; }
ln -sf $A contigs.fa
samtools faidx contigs.fa
step "chromap index + Hi-C mapping"
chromap -i -r contigs.fa -o contigs.index 2> chromap_index.log
chromap --preset hic -r contigs.fa -x contigs.index --remove-pcr-duplicates -1 $R1 -2 $R2 \
  --SAM -o /dev/stdout -t $T 2> chromap.log | samtools view -@ 8 -b - | samtools sort -n -@ 8 -m 2G -T $TMPDIR/sqc.$$ -o hic.bam -
step "yahs"
yahs -o yahs contigs.fa hic.bam > yahs.log 2>&1
rm -f hic.bam   # tens of GB; remap for contact maps if needed
step "gfastats"
for f in contigs.fa yahs_scaffolds_final.fa; do echo "== $f"; gfastats $f 2> /dev/null; done > gfastats.txt
step "compleasm"
compleasm run -a yahs_scaffolds_final.fa -o compleasm -l mammalia --odb odb10 -L $LIN -t $T > compleasm.log 2>&1
cat compleasm/summary.txt
step "merqury"
mkdir -p merqury
cd merqury || exit 1
merqury.sh $MDB ../yahs_scaffolds_final.fa mq > merqury.log 2>&1
cat mq.qv mq.completeness.stats
cd ..
rm -f contigs.index
step done
