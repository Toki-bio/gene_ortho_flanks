"""Build the repeat library for the indel search: RepBase_filtered.bnk (name<TAB>class<TAB>species header, one-line sequence) +
SINEBase.nr95.fa + the five Dip consensuses, as FASTA with RepeatMasker-style >name#class headers."""
import re,sys
BNK,SB,DIP,OUT=sys.argv[1:5]
out=[];seen=set()
def add(n,c,s):
    n=re.sub(r"[^A-Za-z0-9_.\-]","_",n);k=n;i=1
    while k in seen: i+=1;k=f"{n}_{i}"
    seen.add(k);out.append(f">{k}#{c}\n{s}\n")
L=open(BNK,encoding="utf-8",errors="replace").read().split("\n")
for i in range(0,len(L)-1,2):
    if L[i].startswith(">"):
        f=L[i][1:].split("\t");add(f[0],f[1].replace(" ","_") if len(f)>1 else "Unknown",L[i+1].strip().upper())
def fa(p,prefix,cls):
    cur=None;s=[]
    for l in open(p):
        l=l.strip()
        if l.startswith(">"):
            if cur: add(prefix+cur,cls,"".join(s).upper())
            cur=l[1:].split()[0];s=[]
        else: s.append(l)
    add(prefix+cur,cls,"".join(s).upper())
fa(SB,"SB_","SINE/SINEBase"); fa(DIP,"DipSicista_","SINE/Dip")
open(OUT,"w").write("".join(out))
