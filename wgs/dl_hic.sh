#!/bin/bash
# Download ENA Hi-C ERR18010482 (S. trizona, same sample SAMEA120392747 as the HiFi run) into DIR,
# resuming partial files; md5-check; touch DIR/download.ok or DIR/download.failed.
# Usage: dl_hic.sh DIR
set -euo pipefail
H=$1
cd $H || exit 1
[ "$(pwd)" = "$H" ] || exit 1
trap 'touch $H/download.failed' ERR
rm -f download.ok download.failed
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
touch download.ok
echo "download ok $(date '+%F %T')"
