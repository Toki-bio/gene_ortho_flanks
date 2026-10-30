#!/bin/bash
# Download trizona HiFi reads (ENA ERR16751478) in 16 byte ranges, then join + md5 check.
set -euo pipefail
D=/home/toki/gof/reads
mkdir -p $D
cd $D || exit 1
[ "$(pwd)" = "$D" ] || exit 1
U=https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR167/078/ERR16751478/ERR16751478.fastq.gz
SIZE=28552695836
MD5=86da81bc7390f4b24daf22236c5f0c13
N=16
CH=$(( (SIZE + N - 1) / N ))
get_part() {
  local i=$1 a=$(( $1 * CH )) b=$(( ($1 + 1) * CH - 1 ))
  [ $b -ge $SIZE ] && b=$((SIZE - 1))
  local want=$(( b - a + 1 )) f=part_$(printf %02d $i)
  touch $f
  while :; do
    local have; have=$(stat -c %s $f)
    [ "$have" -ge "$want" ] && break
    curl -sL --max-time 3600 -r $((a + have))-$b "$U" >> $f || sleep 5
  done
  echo "part $i done $(stat -c %s $f) / $want"
}
echo "start $(date)"
for i in $(seq 0 $((N - 1))); do get_part $i & done
wait
cat $(ls part_* | sort) > ERR16751478.fastq.gz
got=$(md5sum ERR16751478.fastq.gz | cut -d' ' -f1)
echo "md5 $got (expected $MD5)"
[ "$got" = "$MD5" ] && rm -f part_* && echo "OK $(date)" || echo "MD5 MISMATCH - parts kept"
