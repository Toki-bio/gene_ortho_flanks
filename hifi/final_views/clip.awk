# SAM to small SAM: name gets the phase call, soft clips become hard clips, no quality, no tags. var: calls (file: read call)
BEGIN{ FS=OFS="\t"; while((getline l < calls)>0){ split(l,p," "); C[p[1]]=p[2] } }
/^@/{print; next}
{ c=$6; lead=0; tail=0
  if(match(c,/^[0-9]+S/)){ lead=substr(c,1,RLENGTH-1)+0; c=lead"H" substr(c,RLENGTH+1) }
  if(match(c,/[0-9]+S$/)){ tail=substr(c,RSTART,RLENGTH-1)+0; c=substr(c,1,RSTART-1) tail"H" }
  s=substr($10,lead+1,length($10)-lead-tail)
  k=($1 in C)?C[$1]:"unset"
  print k"|"$1,$2,$3,$4,$5,c,"*",0,0,s,"*" }
