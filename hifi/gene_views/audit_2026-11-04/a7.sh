set -uo pipefail
W=/home/toki/gof/gene_views/TMEM114/audit; cd $W || exit 1
M=/home/toki/miniforge3/envs/hifiasm/bin/minimap2
T1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
ls -la /home/toki/gof/trizona.asm20.mmi | awk '{print $5,$6,$7,$8}'; grep -rl 'trizona.asm20.mmi' /home/toki/gof/*.sh /home/toki/gof/code/*.sh 2>/dev/null | head -3
$M -x asm5 -c --secondary=no -t 24 $T1 l1.fa seg.fa 2>/dev/null > l1t.paf
python3 - <<'PY'
import re
r=[l.rstrip("\n").split("\t") for l in open("/home/toki/gof/geneloss2/final.tsv")]; h=r[0]
st={}
for x in r[1:]:
    d=dict(zip(h,x)); m=re.match(r"(.+):(\d+)-(\d+)",d["hap1_locus"])
    if m: st[f"{m[1]}:{int(m[2])+1}-{m[3]}"]=(d["gene"],d["hap1"].split("(")[0],d["hap1_before_repair"].split("(")[0])
best={}
for l in open("l1t.paf"):
    f=l.split("\t")
    if f[0]=="seg":
        if int(f[10])>=1000: print("seg",f[2],f[3],f[5],f[7],f[8],f[9],f[10])
        continue
    q=f[0]
    if q not in best or int(f[9])>best[q][0]: best[q]=(int(f[9]),int(f[1]),f[5])
import collections
c=collections.Counter()
for q,v in st.items():
    b=best.get(q)
    k="unmapped" if not b else "cov>=0.9" if b[0]/b[1]>=.9 else "cov0.5-0.9" if b[0]/b[1]>=.5 else "cov<0.5"
    same = b and b[2]==q.split(":")[0]
    c[(k,v[2],bool(same))]+=1
for k,v in sorted(c.items(),key=str): print(k,v)
print("unmapped examples",[v[0] for q,v in st.items() if q not in best][:12])
PY
