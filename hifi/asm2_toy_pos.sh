#!/bin/bash
# Positive control for round 2: gaps accepted in round 1 (insert/join/genesearch) that the new sources also span are
# re-planned from the NEW sources only and verified with the round-1 fished reads (which cover those junctions).
# Expect most to be accepted again (different fill source, same reads).
set -euo pipefail
cd /home/toki/gof/asm2/toyrep || exit 1
C=/home/toki/gof/code/hifi
D1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
export PATH=/home/toki/miniforge3/envs/hifiasm/bin:/home/toki/miniforge3/envs/jf/bin:/home/toki/miniforge3/envs/lrtools/bin:$PATH
# round-1 accepted gaps by action, then keep those spanned by hicanu (toy) or wgsp
G=$(python3 - <<'PY'
acc = {}
for l in list(open('/home/toki/gof/repair/dep1/repairs.tsv'))[1:]:
    f = l.rstrip('\n').split('\t')
    if f[13] == 'yes': acc[f[1]] = f[6]
rows = [l.rstrip('\n').split('\t') for l in open('gf_dep1/gaps.tsv')]
ix = {k: i for i, k in enumerate(rows[0])}
pick = {}
for r in rows[1:]:
    a = acc.get(r[0])
    if a and (r[ix['hicanu_status']] == 'spanned' or r[ix['wgsp_status']] == 'spanned'):
        pick.setdefault(a, []).append(r[0])
print(','.join(sum([v[:3] for v in pick.values()], [])))
PY
)
echo "positive-control gaps: $G"
python3 $C/repair_gaps.py plan $D1 gf_dep1/gaps.tsv gf_dep1/fills.fa gf_dep1/fills.fa gf_dep1/fill_cn.tsv /home/toki/gof/gf_all/flankcn/gap_causes.tsv hicanu pos --gaps "$G"
n=$(grep -c '>' pos/cand.fa)
jellyfish count -m 25 -C -s $((4 * n + 1000)) -t 8 --if pos/cand.fa -o pos/cand.jf $D1 pos/pieces.fa
jellyfish dump -c pos/cand.jf > pos/cand.counts
python3 $C/repair_gaps.py words pos pos/cand.counts
python3 $C/repair_gaps.py verify $D1 pos fished.fq.gz dep1 8 --edited-only
rm -f pos/dep1.repaired.fa pos/dep1.genesearch.fa pos/cand.jf pos/cand.counts
cut -f2,6-8,10-15 pos/repairs.tsv
