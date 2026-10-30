#!/bin/bash
cd /home/toki/gof/b2t
export PATH=/home/toki/toga2/env_mlc2/bin:$PATH
Q=/home/toki/toga2/inputs/trizona_nodot.2bit
for g in EPS8L1_237 TMED4_20; do
  grep -P "^$g\t" no_anchor_profile_n.tsv | cut -f11 > /dev/null
  loc=$(python3 /home/toki/gof/code/tests/nb_ncheck.py <(printf "gene\n$g\tx\tx\tx\tx\tx\tx\tx\tx\tx\n") out /home/toki/toga2/b2t/full/toga_out $Q | tail -1 | cut -f11)
  site=$(awk -F'\t' -v g=$g '$1==g{split($7,a,/[:-]/); print a[2]}' out/verdicts.tsv)
  echo "== $g neighbours $loc, empty site at $site"
  python3 /home/toki/gof/code/tests/nrun_pos.py $Q $loc $site
done
