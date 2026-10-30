#!/bin/bash
# Wait for the HiFi download (md5 OK line), then fish/assemble/evaluate the 22 regions.
set -u
W=/home/toki/gof/hifi
cd $W || exit 1
L=/home/toki/gof/reads/download.log
until grep -q "^OK \|MD5 MISMATCH" $L; do sleep 60; done
grep -q "MD5 MISMATCH" $L && { echo "download md5 mismatch - not running"; exit 1; }
for s in fish assemble evaluate; do bash /home/toki/gof/code/hifi/rebuild.sh $s || { echo "stage $s failed"; exit 1; }; done
