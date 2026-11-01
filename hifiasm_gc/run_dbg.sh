#!/bin/bash
# Debug rerun: GC correction ON + HIFIASM_GC_DEBUG on the reads of the R01/R04/R07 piece unitigs (asm25 r_utg)
# and one ROH control unitig. Bins hard-linked from asm_gc (our own copies).
set -euo pipefail
D=/home/toki/gof/hifiasm_gc
O=$D/asm_dbg
mkdir -p $O
cd $O || exit 1
[ "$(pwd)" = "$O" ] || exit 1
export TMPDIR=/home/toki/tmp
for f in tri.ec.bin tri.ovlp.source.bin tri.ovlp.reverse.bin tri.hic.lk.bin tri.hic.tlb.bin; do [ -e $f ] || ln $D/asm_gc/$f $O/$f; done
awk -F'\t' '$1=="A" && ($2=="utg004638l" || $2=="utg005663l" || $2=="utg008327l") {print $5}' /home/toki/gof/asm25/tri.hic.r_utg.noseq.gfa | sort -u > dbg_reads.txt
wc -l dbg_reads.txt
HIFIASM_GC_CORR=on HIFIASM_GC_DEBUG=$O/dbg_reads.txt $D/gc2/hifiasm -o $O/tri -t 48 \
  --h1 /home/toki/gof/hic/ERR18010482_1.fastq.gz --h2 /home/toki/gof/hic/ERR18010482_2.fastq.gz \
  /home/toki/gof/reads/ERR16751478.fastq.gz > hic.log 2>&1
grep -c GCDBG hic.log
for h in hap1 hap2; do cmp -s tri.hic.$h.p_ctg.gfa $D/asm_gc/tri.hic.$h.p_ctg.gfa && echo "$h identical to asm_gc" || echo "$h differs from asm_gc"; done
echo dbg_done
