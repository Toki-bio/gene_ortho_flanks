#!/bin/bash
cd /home/toki/toga2/b2t
ls -la | head -40
echo "== b2t_full.sh"; cat b2t_full.sh
echo "== isoforms"; head -3 sbe_ref_isoforms.tsv; echo "== annot"; head -2 sbe_ref_annotation.bed | cut -f1-6,10
echo "== query_genes.tsv"; head -3 full/toga_out/query_genes.tsv
echo "== orth"; awk -F'\t' '$5=="one2one"' full/toga_out/orthology_classification.tsv | head -2
echo "== loss tx"; grep -m3 "^TRANSCRIPT" full/toga_out/loss_summary.tsv
echo "== tools"; for e in sinederella env_mlc2; do ls -d /home/toki/*/envs/$e /home/toki/miniforge3/envs/$e /home/toki/miniconda3/envs/$e 2>/dev/null; done
which minimap2 samtools bedtools python3 2>&1; conda env list 2>/dev/null | head -20
df -h /home | tail -1
