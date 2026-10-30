#!/bin/bash
cd /home/toki/gof/toy || exit 1
ls
for v in keep del gap; do echo "== $v: gene call"; [ -s out_$v/verdicts.tsv ] && awk -F'\t' 'NR>1{print $NF}' out_$v/verdicts.tsv | sort | uniq -c; done
ls -la /home/toki/gof/reads | head -5
