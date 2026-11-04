import sys,re,collections
m8,chunks=sys.argv[1],sys.argv[2]
L={};n=None
for l in open("C:/Temp/gv/lib.fa"):
    if l.startswith(">"): n=l[1:].strip(); L[n]=0
    else: L[n]+=len(l.strip())
qlen={};n=None
for l in open(chunks):
    if l.startswith(">"): n=l[1:].strip(); qlen[n]=0
    else: qlen[n]+=len(l.strip())
hits=collections.defaultdict(list)
for l in open(m8):
    f=l.rstrip("\r\n").split("\t")
    if len(f)<12 or l.startswith("#"): continue
    q,s=f[0],f[1]; idn=float(f[2]); qs,qe=sorted((int(f[6]),int(f[7]))); ss,se=sorted((int(f[8]),int(f[9]))); bits=float(f[11])
    hits[q].append((bits,qs,qe,s,idn,ss,se))
def base(q): return q.split("@")[0], int(q.split("@")[1]) if "@" in q else 0
elem=collections.defaultdict(list)
for q,h in hits.items():
    h.sort(reverse=True); acc=[]
    for b,qs,qe,s,idn,ss,se in h:
        if qe-qs+1<80 or b<50: continue
        if any(min(qe,a[2])-max(qs,a[1])>0.3*(qe-qs+1) for a in acc): continue
        acc.append((b,qs,qe,s,idn,ss,se))
    nm,off=base(q)
    for b,qs,qe,s,idn,ss,se in acc:
        cls=s.split("#")[1].split("/")[0] if "#" in s else "?"
        full=s
        elem[nm].append((qs+off,qe+off,s,idn,ss,se,L.get(s,0),b))
tot=collections.Counter(); 
allq={}
for q,l in qlen.items():
    nm,off=base(q); allq[nm]=allq.get(nm,0)+l
for nm in sorted(elem):
    for qs,qe,s,idn,ss,se,sl,b in sorted(elem[nm]):
        cls=s.split("#")[1]; 
        print(f"{nm}\t{qs}-{qe}\t{qe-qs+1}\t{s}\t{idn:.0f}%\t{ss}-{se}/{sl}\t{b:.0f}")
