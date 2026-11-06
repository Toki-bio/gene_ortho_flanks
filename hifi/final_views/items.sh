#!/bin/bash
# what gets a view: all new fills, every gene gap still open, the gene that changed, the median gap of each cause per haplotype, TMED4 as control
set -euo pipefail
. "$(dirname "$0")/cfg.sh"
cd "$W" || exit 1
S=$F/gap_status.tsv
# genes per gap: hap1 gene-search still gap, or changed since geneloss2
awk -F'\t' -v OFS='\t' 'FILENAME==ARGV[1]{ p[$1]=$15; next } FNR==1{next}
	{ ch=(p[$1]~/gap/ && $15!~/gap/); if($15~/gap/||ch){ split($6,a,/[(:]/); print a[2],$1 } }' $G/geneloss2/final.tsv $GL/final.tsv \
	| sort | awk -F'\t' '{ if($1!=g){ if(g!="")print g"\t"l; g=$1; l=$2 } else l=l","$2 } END{print g"\t"l}' > gene_gaps.tsv
awk -F'\t' -v OFS='\t' '
	FILENAME==ARGV[1]{ GG[$1]=$2; next }
	FNR==1{next}
	{ id=$1"_"$2; k=""; lab=$11
	  if($7=="lr") k="new"
	  if($1=="dep1" && ($2 in GG)){ k=k (k==""?"":"+") "gene"; lab=GG[$2] }
	  if($1=="dep1" && $2=="gap0441"){ k=k (k==""?"":"+") "control"; lab="TMED4_20" }
	  if(k!="") print id,k,$1,$2,$3,$4,$5,lab,$6,$8,$9,$10,$11
	  else if($6=="N"){ c=$1 SUBSEP $11; n[c]++; L[c,n[c]]=id OFS "blocked" OFS $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $11 OFS $6 OFS $8 OFS $9 OFS $10 OFS $11 } }
	END{ for(c in n) print L[c,int((n[c]+1)/2)] }' gene_gaps.tsv "$S" | sort > items.tsv
cut -f2 items.tsv | sort | uniq -c
