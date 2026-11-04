#!/bin/bash
# QC of the HiCanu contigs: compleasm (mammalia odb10, offline library) and Merqury QV/completeness against the HiFi
# k-mer db, same calls as wgs/run_sqc.sh. The two steps do not depend on each other: a failed compleasm must not
# skip Merqury. HiCanu holds both parental copies, so read duplication in compleasm as expected, not as a defect.
# Usage: asm2_canu_qc.sh OUT_DIR CONTIGS.fa [THREADS]
set -uo pipefail
O=$1
A=$2
T=${3:-16}
E=/home/toki/miniforge3/envs/asm25
export PATH=$E/bin:$PATH
export TMPDIR=/home/toki/tmp
export MERQURY=$(ls -d $E/share/merqury*)
mkdir -p "$O/merqury"
cd "$O" || exit 1
echo "#### compleasm $(date '+%F %T')"
compleasm run -a "$A" -o compleasm -l mammalia --odb odb10 -L /home/toki/gof/compleasm_lib -t "$T" > compleasm.log 2>&1
echo "compleasm exit $?"
cat compleasm/summary.txt
echo "#### merqury $(date '+%F %T')"
cd merqury || exit 1
merqury.sh /home/toki/gof/sqc/hifi.meryl "$A" mq > merqury.log 2>&1
echo "merqury exit $?"
cat mq.qv mq.completeness.stats
echo "#### done $(date '+%F %T')"
