# local: alignment structure
def rd(p):
    d={};n=None
    for l in open(p):
        l=l.strip()
        if l.startswith(">"): n=l[1:].split()[0]; d[n]=[]
        else: d[n].append(l)
    return {k:''.join(v).upper() for k,v in d.items()}
a=rd("C:/work/gene_ortho_flanks/views/het_genes/TMEM114.aln.fa")
n=len(next(iter(a.values())))
for k,v in a.items():
    c=[i for i,x in enumerate(v) if x!='-']; print(k,"first",c[0]+1,"last",c[-1]+1,"bases",len(c),"N",v.count("N"))
ks=list(a)
def ident(x,y,s,e):
    m=t=0
    for i in range(s,e):
        p,q=a[x][i],a[y][i]
        if p!='-' and q!='-': t+=1; m+=p==q
    return (m/t if t else 0,t)
print("window identity tri1/tri2 and tri2/bet1")
for s in range(0,n,5000):
    e=min(n,s+5000); i1=ident(ks[0],ks[1],s,e); i2=ident(ks[1],ks[2],s,e); i3=ident(ks[2],ks[3],s,e)
    print(s+1,e,f"t1t2 {i1[0]:.3f}/{i1[1]}  t2b1 {i2[0]:.3f}/{i2[1]}  b1b2 {i3[0]:.3f}/{i3[1]}")
ex=[l.split("\t") for l in open("C:/work/gene_ortho_flanks/views/het_genes/TMEM114.bed")]
for f in ex:
    if f[0]=="TMEM114_bet_hap1" and "exon" in f[3]:
        s,e=int(f[1]),int(f[2]); print(f[3],"t2/b1",ident(ks[1],ks[2],s,e),"hap1 bases",sum(c!='-' for c in a[ks[0]][s:e]))
print(sorted(set(l.split("\t")[0] for l in open("C:/work/gene_ortho_flanks/views/het_genes/TMEM114.bed"))))
