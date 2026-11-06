#!/bin/bash
# all stages. Argument toy = three items first
set -uo pipefail
cd /home/toki/gof/final_views/code || exit 1
bash items.sh || exit 1
cd .. || exit 1
if [ "${1:-}" = toy ]; then
	grep -P '^dep1_gap0441\t|^dep1_gap0494\t|^dep1_gap0086\t' items.tsv > toy.tsv
	L=toy.tsv
else L=items.tsv; fi
rm -rf w out
bash code/win.sh $L || exit 1
bash code/fish.sh || exit 1
cut -f1 $L | xargs -P 12 -L 1 bash code/one.sh
echo finished
