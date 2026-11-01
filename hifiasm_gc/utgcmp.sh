#!/bin/bash
cd /home/toki/gof/hifiasm_gc/eval || exit 1
paste utg_rule_asm25.tsv utg_rule_gc.tsv | awk -F'\t' 'NR==1 || ($1<=24) {print $1"\t"$2"\t"$3"\t"$4"\tasm25_both="$5"\tgc_both="$12"\tgc_one="$13}'
