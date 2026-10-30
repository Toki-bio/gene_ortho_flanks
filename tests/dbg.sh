#!/bin/bash
ps -u toki -o pid,ppid,etime,args | grep -E "curl|minimap2|gof.py|run_toy|get_trizona" | grep -v grep | cut -c1-200
cd /home/toki/tmp || exit 1
curl -sL -r 1784543490-1784543999 -o t.bin -w "%{http_code} %{size_download}\n" --max-time 30 https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR167/078/ERR16751478/ERR16751478.fastq.gz; rm -f t.bin
