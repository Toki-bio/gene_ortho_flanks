#!/bin/bash
cd /home/toki/tmp || exit 1
u=https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR167/078/ERR16751478/ERR16751478.fastq.gz
start=$(date +%s.%N)
for i in 0 1 2 3 4 5 6 7; do
  a=$((i * 1000000000)); b=$((a + 299999999))
  curl -sL -r $a-$b -o /dev/null -w "%{size_download}\n" --max-time 40 "$u" > sp_$i.txt &
done
wait
end=$(date +%s.%N)
tot=$(cat sp_*.txt | awk '{s+=$1} END{print s}')
python3 -c "print('8 parallel: %.1f MB in %.0f s = %.2f MB/s' % ($tot/1e6, $end-$start, $tot/1e6/($end-$start)))"
rm -f sp_*.txt
