#!/bin/bash
# one row per item from out/, and the fresh verdict against the one in the final table
set -uo pipefail
. "$(dirname "$(readlink -f "$0")")/cfg.sh"
cd "$W" || exit 1
printf 'id\tkind\tasm\tgap\tcontig\tstart1\tend1\tlabel\tstate\taction\tverdict_final\tlr_status\tlr_cause\treads\tT\tO\ttie\tforeign\tshort\tone_side\tT_cross\tO_cross\ttie_cross\tT_span\tT_x\tT_x_median\tT_indel_median\tT_big\tO_x_median\ttie_x_median\tT_enter\tT_enter_x\tO_enter_x\tverdict_fresh\tfill_len\tview_reads\tsupport\tbet\tother\n' > views.tsv
while IFS=$'\t' read -r id kind a gap c s e lab st act v lrs lrc; do
	o=out/$id
	if [ ! -s $o/sum.tsv ]; then echo "no output: $id"; continue; fi
	. w/$id/meta
	sup=$(awk -F'\t' -v a=$a -v g=$gap '$1==a&&$6==g{print $12}' $LR/new_fills.tsv)
	bet=$(cat $o/bet.loc 2>/dev/null | tr ' ' ':')
	oth=$(tr '\n' ';' < w/$id/oreg 2>/dev/null)
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' "$id" "$kind" "$a" "$gap" "$c" "$s" "$e" "$lab" "$st" "$act" "$v" "$lrs" "$lrc" >> views.tsv
	tr -d '\n' < $o/sum.tsv >> views.tsv
	printf '\t%s\t%s\t%s\t%s\t%s\n' $(( v1 - v0 )) "$(samtools view -c $o/view.bam)" "$sup" "$bet" "$oth" >> views.tsv
done < <(sed 's/		/	.	/g;s/		/	.	/g' items.tsv)
echo "rows $(( $(wc -l < views.tsv) - 1 ))"
echo "== new fills, final verdict against fresh reads"
awk -F'\t' 'NR>1 && $2~/new/{print $11" -> "$34}' views.tsv | sort | uniq -c
du -sh out; find out -name aln.fa -size 0 | head; ls -la views.tsv
