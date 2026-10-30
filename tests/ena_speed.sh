#!/bin/bash
cd /home/toki/tmp || exit 1
for u in https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR167/078/ERR16751478/ERR16751478.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR167/078/ERR16751478/ERR16751478.fastq.gz; do
  curl -sL -r 0-299999999 -o /dev/null -w "%{http_code} bytes %{size_download} in %{time_total} s = %{speed_download} B/s\n" --max-time 60 "$u"
done
ls -d /home/toki/miniforge3/envs/*/bin/hifiasm /home/toki/miniforge3/envs/*/bin/flye /home/toki/miniforge3/bin/hifiasm 2>/dev/null
which hifiasm flye 2>/dev/null
