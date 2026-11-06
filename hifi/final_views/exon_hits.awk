# ssearch36 -m 8 to BED: best hit per exon and row, at least 70 % of the exon, E below 1e-3. Query names GENE_exonN_LEN
BEGIN{ FS=OFS="\t" }
NF>=12 && !/^#/{ split($1,p,"_"); L=p[length(p)]+0; a=$9<$10?$9:$10; b=$9<$10?$10:$9; n=$8-$7; if(n<0)n=-n; n++; k=$1 SUBSEP $2
	if(!(k in B)||$12>B[k]){B[k]=$12; A[k]=a; Z[k]=b; N[k]=n; E[k]=$11; I[k]=$3; LL[k]=L} }
END{ for(k in B) if(N[k]>=0.7*LL[k] && E[k]<1e-3){ split(k,q,SUBSEP); m=q[1]; sub(/_[0-9]+$/,"",m); sub(/_exon/," exon",m)
	print q[2],A[k]-1,Z[k],m" ("LL[k]" bp, "int(I[k]+0.5)"%)",0,"+",A[k]-1,Z[k],"0,0,150" } }
