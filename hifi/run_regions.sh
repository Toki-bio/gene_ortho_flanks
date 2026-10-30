#!/bin/bash
set -euo pipefail
cd /home/toki/gof/b2t || exit 1
export PATH=/home/toki/toga2/env_mlc2/bin:$PATH
mkdir -p /home/toki/gof/hifi
python3 /home/toki/gof/code/hifi/make_regions.py out no_anchor_profile_n.tsv /home/toki/toga2/inputs/trizona_nodot.2bit 3 > /home/toki/gof/hifi/regions.tsv
column -t /home/toki/gof/hifi/regions.tsv
