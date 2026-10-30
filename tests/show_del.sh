#!/bin/bash
cd /home/toki/gof/b2t
awk -F"	" 'NR==1 || $10=="deletion_larger_than_flanks"' no_anchor_profile_n.tsv | cut -f1,4,5,9,11,12,13 | column -t
