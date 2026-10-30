#!/bin/bash
# Toy for rebuild.sh: simulated HiFi reads with a known answer (see make_sim_toy.py).
set -euo pipefail
T=/home/toki/gof/hifi_toy
mkdir -p $T
cd $T || exit 1
export PATH=/home/toki/toga2/env_mlc2/bin:$PATH
python3 /home/toki/gof/code/hifi/make_sim_toy.py /home/toki/gof/hifi/regions.tsv /home/toki/toga2/inputs/trizona_nodot.2bit \
  /home/toki/toga2/b2t/sbe.2bit $T R00:insert R06:empty R10:asis R19:control R20:control R21:control
W=$T READS=$T/reads.fastq.gz bash /home/toki/gof/code/hifi/rebuild.sh all 2>&1 | grep -v "^\[M::"
echo "== truth"; for f in truth_*.txt; do echo "$f $(cat $f)"; done
echo "background reads fished: $(grep -c '^@bg_' fished.tsv || true)"
