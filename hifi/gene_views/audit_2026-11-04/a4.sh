cd /home/toki/gof || exit 1
W=/home/toki/gof/gene_views/TMEM114/audit
E=/home/toki/miniforge3/envs; S=$E/sinederella/bin/samtools; M=$E/hifiasm/bin/minimap2
B1=/home/toki/Sicista2026/assembly/hifiasm_out/sicista_hifiasm.p_ctg.fa
echo "== recall header"; sed -n 4,13p geneloss2/recall.sh | cut -c1-170
echo "== prep"; grep -n -E '2bit|\.fa|bed' geneloss2/run.sh | cut -c1-170 | head -12; ls geneloss2/prep | head
echo "== B1 names"; head -3 $B1.fai | cut -f1,2; wc -l $B1.fai
echo "== files holding ptg000020l"
for f in /home/toki/Sicista2026/assembly/*/*.fa /home/toki/Sicista2026/assembly/*.fa /home/toki/Sicista2026/reference/*.fa /home/toki/Sicista2026/reference/*.fna; do [ -f $f.fai ] && grep -l -P '^ptg000020l\t' $f.fai; done 2>/dev/null | head
echo "== bet rows vs B1"
$M -x asm10 -c --secondary=no -t 16 $B1 $W/rows.fa 2>/dev/null | awk '$10>2000{print $1,$5,$3,$4,$6,$7,$8,$9,$10,$11}'
echo "== pairwise rows"
cd $W
python3 - <<'PY'
import subprocess,re
M="/home/toki/miniforge3/envs/hifiasm/bin/minimap2"
seq={};n=None
for l in open("rows.fa"):
    l=l.strip()
    if l[0]=='>': n=l[1:]; seq[n]=''
    else: seq[n]+=l
for k,v in seq.items(): open(k+".fa","w").write(f">{k}\n{v}\n")
ks=list(seq)
pairs=[(ks[0],ks[1]),(ks[1],ks[2]),(ks[1],ks[3]),(ks[2],ks[3])]
for t,q in pairs:
    o=subprocess.run([M,"-x","asm20","-c","--cs=short","-r","2000,120000","-z","1000000","--secondary=no",t+".fa",q+".fa"],capture_output=True,text=True).stdout
    print("#",t,q)
    for l in o.splitlines():
        f=l.split("\t"); print(" block q",f[2],f[3],f[4],"t",f[7],f[8],"m",f[9],f[10])
        cs=[x for x in f if x.startswith("cs:Z:")][0][5:]
        qp=int(f[2]); tp=int(f[7])
        if f[4]=='-': print("  minus strand, skip"); continue
        for op,arg in re.findall(r"([:*+\-~])([A-Za-z0-9]+)",cs):
            if op==':': qp+=int(arg); tp+=int(arg)
            elif op=='*': qp+=1; tp+=1
            elif op=='+':
                if len(arg)>=50: print(f"  INS_in_{q} q{qp+1}-{qp+len(arg)} at t{tp} len {len(arg)}")
                qp+=len(arg)
            elif op=='-':
                if len(arg)>=50: print(f"  INS_in_{t} t{tp+1}-{tp+len(arg)} at q{qp} len {len(arg)}")
                tp+=len(arg)
PY
echo "== N-gene candidates"
cd /home/toki/gof; grep -E '^(EPS8L1|TMED4)' geneloss2/final.tsv | cut -f1,3,6,7,14,15,16 | head
grep -c . gf_all/flankcn/gap_genes.tsv; grep -w -E 'EPS8L1|TMED4' gf_all/flankcn/gap_genes.tsv | cut -c1-140 | head -4
