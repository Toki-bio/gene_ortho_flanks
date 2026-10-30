#!/bin/bash
set -euo pipefail
cd /home/toki/gof/b2t || exit 1
export PATH=/home/toki/toga2/env_mlc2/bin:/home/toki/miniforge3/envs/sinederella/bin:$PATH
T=/home/toki/toga2/b2t
python3 /home/toki/gof/code/make_aln.py out $T/sbe_ref_annotation.bed $T/sbe_ref_isoforms.tsv $T/sbe.2bit /home/toki/toga2/inputs/trizona_nodot.2bit EPS8L1_237 TMED4_20
