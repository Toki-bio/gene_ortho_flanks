import re,csv,collections
def rd(p):
    d={};n=None
    for l in open(p):
        l=l.strip()
        if l.startswith(">"): n=l[1:].split()[0]; d[n]=[]
        else: d[n].append(l)
    return {k:''.join(v).upper() for k,v in d.items()}
a=rd("C:/work/gene_ortho_flanks/views/het_genes/TMEM114.aln.fa")
BASE={"TMEM114_tri_hap1":("OZ418350.1",111012847),"TMEM114_tri_hap2":("OZ418090.1",110274646)}
def runs(x,y,gap):
    rx,ry=a[x],a[y]; n=len(rx); R=[]; i=0
    while i<n:
        s=1 if (rx[i]!='-' and ry[i]=='-') else 2 if (ry[i]!='-' and rx[i]=='-') else 0
        if s:
            j=i
            while j<n and ((rx[j]!='-' and ry[j]=='-') if s==1 else (ry[j]!='-' and rx[j]=='-')): j+=1
            R.append((s,i,j)); i=j
        else: i+=1
    m=[]
    for s,i,j in R:
        if m and m[-1][0]==s and i-m[-1][2]<=gap: m[-1]=(s,m[-1][1],j)
        else: m.append((s,i,j))
    out=[]
    for s,i,j in m:
        rn=x if s==1 else y; seq=a[rn][i:j].replace("-","")
        if len(seq)>=50: out.append((rn,i,j,len(a[rn][:i].replace("-",""))))
    return out
def coord(rn,idx0,idx1):
    if rn in BASE:
        c,b=BASE[rn]; return f"{c}:{b+idx0:,}-{b+idx1-1:,}"
    return f"{rn.split('_')[-2]}_{rn.split('_')[-1]} row {idx0+1}-{idx1}"
tables=[]
for tag,(x,y,gap,fn) in {"hap1_vs_hap2":("TMEM114_tri_hap1","TMEM114_tri_hap2",20,"cls_tri12.dedup.tsv"),"trizona_vs_betulina":("TMEM114_tri_hap2","TMEM114_bet_hap1",100,"cls_TB.dedup.tsv")}.items():
    ind=runs(x,y,gap)
    # names as used earlier
    if tag=="hap1_vs_hap2":
        pieces=[]; names={}
        k=0
        for (rn,i,j,u) in ind:
            k+=1; names[f"I{k:03d}"]=(rn,i,j,u)
        big=[names[f"I{n:03d}"] for n in range(1,17)]
        # element offset within concatenated big insertion -> piece
        cum=[];o=0
        for rn,i,j,u in big:
            L=len(a[rn][i:j].replace("-","")); cum.append((o,o+L,rn,u)); o+=L
        look={"INS65kb_hap2only_I001-I016":("big",cum)}
        for n in range(17,len(ind)+1):
            rn,i,j,u=names[f"I{n:03d}"]; look[[k2 for k2 in []] and "" or f"I{n:03d}"]=("one",(rn,u))
    else:
        look={}
        for k,(rn,i,j,u) in enumerate(ind,1): look[f"TB{k:03d}"]=("one",(rn,u))
    rows=list(csv.reader(open(fn),delimiter="\t"))[1:]
    for r in rows:
        nm,s,e=r[0],int(r[1]),int(r[2]); key=nm if nm in look else nm.split("_")[0]
        kind,info=look[key]
        if kind=="one": rn,u=info; c=coord(rn,u+s-1,u+e)
        else:
            c="?"
            for o0,o1,rn,u in info:
                if o0<s<=o1: c=coord(rn,u+(s-1-o0),u+(min(e,o1)-o0)); break
        tables.append((tag,nm,c,r[3],r[4],r[5],r[6],r[7]))
with open("TMEM114_indel_TEs.tsv","w") as o:
    o.write("comparison\tindel\tlocation\tlen\tbest_match\tclass\t%id\tmatch_range_in_library_seq\n")
    for t in tables: o.write("\t".join(t)+"\n")
print(len(tables))
for t in tables[:4]+tables[-4:]: print(t)
