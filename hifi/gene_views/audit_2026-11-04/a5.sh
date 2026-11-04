set -uo pipefail
W=/home/toki/gof/gene_views/TMEM114/audit; cd $W || exit 1
E=/home/toki/miniforge3/envs; S=$E/sinederella/bin/samtools; M=$E/hifiasm/bin/minimap2; SS=$E/sinederella/bin/ssearch36
T1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
T2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
tr -d '\r' < lib.crlf.fa > lib.fa
python3 - <<'PY'
seq={};n=None
for l in open("rows.fa"):
    l=l.strip()
    if l[0]=='>': n=l[1:].replace("TMEM114_",""); seq[n]=''
    else: seq[n]+=l
Q=[("tri_hap2",15709,79274,"hap1_deletion"),
("bet_hap1",1593,1757,"tri_vs_bet"),("tri_hap2",3193,3381,"tri_vs_bet"),("bet_hap1",10820,11494,"tri_vs_bet"),("bet_hap1",17667,17902,"bet1_only"),
("bet_hap1",24328,25745,"tri_vs_bet"),("tri_hap2",39941,40248,"tri_vs_bet"),("bet_hap1",44879,44955,"tri_vs_bet"),("tri_hap2",57068,57123,"tri_vs_bet"),
("tri_hap2",58303,59088,"tri_vs_bet"),("tri_hap2",61489,61612,"tri_vs_bet"),("bet_hap1",63245,63435,"tri_vs_bet"),("tri_hap2",68317,68686,"tri_vs_bet"),
("bet_hap1",80517,80781,"tri_vs_bet"),("bet_hap1",81532,86596,"bet1_only"),("tri_hap2",83210,83719,"tri_vs_bet"),("bet_hap1",90664,90902,"tri_vs_bet"),
("bet_hap1",92072,92603,"tri_vs_bet"),("bet_hap2",30905,30987,"tri_vs_bet"),("tri_hap2",32682,32732,"tri_vs_bet"),("bet_hap2",53167,53868,"bet2_only"),
("tri_hap2",83783,83909,"tri_vs_bet"),("bet_hap1",89371,89454,"bet1_only"),("bet_hap2",89783,90091,"bet2_only")]
with open("q.fa","w") as o:
    for r,s,e,c in Q:
        x=seq[r][s-1:e]
        for off in range(0,len(x),5000):
            p=x[off:off+5300]
            if len(p)>=50: o.write(f">{c}|{r}|{s}|{e}|{off}\n{p}\n")
        if len(x)>5000: pass
open("seg.fa","w").write(">seg\n"+seq["tri_hap2"][15708:79274]+"\n")
PY
grep -c '>' q.fa
$SS -n -m 8 -E 1e-5 -T 16 q.fa lib.fa > q.m8 2> q.err; echo ssearch $? $(wc -l < q.m8)
echo "== segment elsewhere on hap1 (blocks >=1 kb)"
$M -x asm20 -c -N 20 -p 0.3 -t 16 /home/toki/gof/trizona.asm20.mmi seg.fa 2>/dev/null | awk '$11>=1000{print $3,$4,$5,$6,$8,$9,$10,$11}' | sort -n | head -30
echo "== tail reads"
python3 - <<'PY'
import subprocess,re
S="/home/toki/miniforge3/envs/sinederella/bin/samtools"
A=[]
for l in subprocess.run([S,"view","../reads.bam"],capture_output=True,text=True).stdout.splitlines():
    f=l.split("\t"); p=int(f[3]); sp=sum(int(n) for n,o in re.findall(r"(\d+)([MDN=X])",f[5])); A.append((f[2],p,p+sp,int(f[1])))
for ref,pos in (("tri_hap2",120246),("tri_hap1",56616)):
    print(ref,pos,"spanning+-1kb",sum(1 for a in A if a[0].startswith(ref) and a[1]<=pos-1000 and a[2]>=pos+1000 and not a[3]&0x800))
PY
echo "== all final.tsv loci back on deposited"
cd /home/toki/gof
python3 - <<'PY'
import re
r=[l.rstrip("\n").split("\t") for l in open("geneloss2/final.tsv")]; h=r[0]
o={1:open("gene_views/TMEM114/audit/l1.bed","w"),2:open("gene_views/TMEM114/audit/l2.bed","w")}
n=0
for x in r[1:]:
    d=dict(zip(h,x))
    for k in (1,2):
        m=re.match(r"(.+):(\d+)-(\d+)",d[f"hap{k}_locus"])
        if m: o[k].write(f"{m[1]}:{int(m[2])+1}-{m[3]}\n"); n+=1
print("loci",n)
PY
cd $W
$S faidx -r l1.bed /home/toki/gof/repair/dep1.repaired.fa > l1.fa; $S faidx -r l2.bed /home/toki/gof/repair/dep2.repaired.fa > l2.fa
$M -x asm5 -c --secondary=no -t 16 /home/toki/gof/trizona.asm20.mmi l1.fa 2>/dev/null > l1.paf
$M -x asm5 -c --secondary=no -t 24 $T2 l2.fa 2>/dev/null > l2.paf
python3 - <<'PY'
for k in (1,2):
    L={}
    for l in open(f"l{k}.bed"): L[l.strip()]=0
    best={}
    for l in open(f"l{k}.paf"):
        f=l.split("\t"); q=f[0]; 
        if q not in best or int(f[9])>best[q][0]: best[q]=(int(f[9]),int(f[1]),f[5],int(f[10]))
    same=sum(1 for q,b in best.items() if b[2]==q.split(":")[0])
    cov=[b[0]/b[1] for b in best.values()]
    print("hap",k,"loci",len(L),"mapped",len(best),"same chromosome",same,"cov>=0.9",sum(c>=.9 for c in cov),"cov<0.5",sum(c<.5 for c in cov))
PY
