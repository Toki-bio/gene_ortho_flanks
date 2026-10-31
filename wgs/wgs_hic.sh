#!/bin/bash
# Hi-C phased hifiasm of S. trizona: download ENA Hi-C ERR18010482 (same sample SAMEA120392747 as the
# HiFi run), then rerun hifiasm with the same -o prefix so it reuses tri.ec.bin / tri.ovlp.*.bin from
# wgs_hifiasm.sh and only does the Hi-C phasing. Output: tri.hic.hap1.p_ctg.fa, tri.hic.hap2.p_ctg.fa.
set -euo pipefail
W=/home/toki/gof/wgs
H=/home/toki/gof/hic
mkdir -p $H
cd $H || exit 1
[ "$(pwd)" = "$H" ] || exit 1
export PATH=/home/toki/miniforge3/envs/hifiasm/bin:$PATH
export TMPDIR=/home/toki/tmp
U=https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR180/082/ERR18010482
echo "download start $(date '+%F %T')"
wget -q -c $U/ERR18010482_1.fastq.gz &
P1=$!
wget -q -c $U/ERR18010482_2.fastq.gz &
P2=$!
wait $P1
wait $P2
printf '%s  %s\n' b708c426184f8f90c20bdbdbca9473ef ERR18010482_1.fastq.gz \
  94dcce3a3465e8714bb22b81997e7ae3 ERR18010482_2.fastq.gz | md5sum -c
echo "download ok $(date '+%F %T')"
cd $W || exit 1
[ -s tri.ec.bin ] && [ -s tri.ovlp.source.bin ] && [ -s tri.ovlp.reverse.bin ] || { echo "bins missing"; exit 1; }
echo "hifiasm hic start $(date '+%F %T') $(hifiasm --version)"
hifiasm -o $W/tri -t 64 -l0 --h1 $H/ERR18010482_1.fastq.gz --h2 $H/ERR18010482_2.fastq.gz \
  /home/toki/gof/reads/ERR16751478.fastq.gz > $W/hifiasm_hic.log 2>&1
for p in hic.hap1.p_ctg hic.hap2.p_ctg; do
  awk '$1=="S"{print ">"$2"\n"$3}' $W/tri.$p.gfa > $W/tri.$p.fa
done
echo "done $(date '+%F %T')"
