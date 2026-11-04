import collections,re
sp={}
L=open("C:/work/MSA-viewer-new/RepBase_filtered.bnk",encoding="utf-8",errors="replace")
for l in L:
    if l.startswith(">"):
        f=l[1:].rstrip("\n").split("\t"); sp[re.sub(r"[^A-Za-z0-9_.\-]","_",f[0])]=f[2] if len(f)>2 else ""
ln={};n=None
for l in open("lib.fa"):
    if l[0]=='>': n=l[1:].strip(); ln[n]=0
    else: ln[n]+=len(l.strip())
H=collections.defaultdict(list)
for l in open("q.m8"):
    f=l.rstrip("\r\n").split("\t")
    if len(f)<12 or l[0]=='#': continue
    c,r,s,e,off=f[0].split("|"); s,e,off=int(s),int(e),int(off)
    a,b=sorted((int(f[6]),int(f[7]))); x,y=sorted((int(f[8]),int(f[9])))
    H[(c,r,s,e)].append((float(f[11]),off+a,off+b,f[1],float(f[2]),x,y,float(f[10])))
BASE={"tri_hap2":lambda p:f"OZ418090.1:{110274645+p-1:,}","bet_hap1":lambda p:f"h1tg000312l:{6327124+p:,}","bet_hap2":lambda p:f"h2tg001095l:{722000-p+1:,}(-)"}
out=[]
for k,h in H.items():
    h.sort(reverse=True); acc=[]
    for t in h:
        if t[2]-t[1]+1<80 or t[0]<50: continue
        if any(min(t[2],a[2])-max(t[1],a[1])>0.3*(t[2]-t[1]+1) for a in acc): continue
        acc.append(t)
    for b,a1,a2,s,idn,x,y,ev in sorted(acc,key=lambda t:t[1]):
        nm,cl=s.split("#"); spc=sp.get(nm,"SINEBase/Dip")
        p=k[2]+a1-1; out.append((k[0],k[1],f"{k[2]}-{k[3]}",k[3]-k[2]+1,BASE[k[1]](p),a2-a1+1,nm,cl,spc,idn,f"{x}-{y}/{ln.get(s,0)}",b,ev))
with open("TMEM114_indel_TEs.tsv","w",newline="\n") as o:
    o.write("comparison\tsequence\tindel_row_coords\tindel_len\tstart_position\thit_len\tbest_match\tclass\tmatch_species\tpct_id\tmatch_range\tbits\tevalue\n")
    for t in out: o.write("\t".join(map(str,t))+"\n")
for c in sorted(set(t[0] for t in out)):
    cnt=collections.Counter();bp=collections.Counter()
    for t in out:
        if t[0]==c:
            g=t[7]; cnt[g]+=1; bp[g]+=t[5]
    print(c,sum(bp.values()),[(g,cnt[g],bp[g]) for g,_ in bp.most_common()])
print("-- non-mammal or weak")
for t in out:
    if t[8] and not re.search("Mus|Rattus|Homo|Mammal|Rodent|Canis|Bos|Sus|Equus|Cavia|Oryct|Eutheria|Primat|Sciur|Jaculus|Cricet|Peromyscus|Sicista|Dip|SINEBase|Loxodonta|Felis",t[8]): print(t[0],t[4],t[5],t[6],t[7],t[8],t[9],t[12])
print("-- non-deletion rows")
for t in out:
    if t[0]!="hap1_deletion": print(t[0],t[1],t[2],"len",t[3],"|",t[5],t[6],t[7],t[8],t[9],t[10])
print("-- deletion >=400"); 
for t in out:
    if t[0]=="hap1_deletion" and t[5]>=400: print(t[4],t[5],t[6],t[7],t[8],t[9],t[10])
