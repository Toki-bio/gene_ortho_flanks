cd /home/toki/gof || exit 1
W=/home/toki/gof/gene_views/TMEM114/audit; mkdir -p $W
E=/home/toki/miniforge3/envs; S=$E/sinederella/bin/samtools; M=$E/hifiasm/bin/minimap2
T1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
T2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
B1=/home/toki/Sicista2026/assembly/hifiasm_out/sicista_hifiasm.p_ctg.fa
B2=/home/toki/gof/betulina/gw/bet_hap2.fa
R1=repair/dep1.repaired.fa; R2=repair/dep2.repaired.fa
tr -d '\r' < gene_views/TMEM114/indels/rows.fa > $W/rows.fa
echo "== what recall used"; cut -c1-170 geneloss2/recall.sh | grep -n -E 'repair|\.fa|dep[12]' | head -12
for f in $R1 $R2 $T1 $T2 $B1 $B2; do [ -f $f.fai ] || $S faidx $f; done
echo "== A rows on repaired chr10"
$S faidx $R1 OZ418350.1 > $W/r1c10.fa; $S faidx $R2 OZ418090.1 > $W/r2c10.fa
$M -x asm10 -c --secondary=no $W/r1c10.fa $W/rows.fa 2>/dev/null | awk '$1~/tri_hap1/&&$10>500{print $1,$5,$3,$4,$6,$8,$9,$10,$11}'
$M -x asm10 -c --secondary=no $W/r2c10.fa $W/rows.fa 2>/dev/null | awk '$1~/tri_hap2/&&$10>500{print $1,$5,$3,$4,$6,$8,$9,$10,$11}'
echo "== B bet rows on betulina assemblies"
$S faidx $B1 ptg000020l > $W/b1c.fa
$M -x asm10 -c --secondary=no $W/b1c.fa $W/rows.fa 2>/dev/null | awk '$10>2000{print $1,$5,$3,$4,$6,$8,$9,$10,$11}'
$M -x asm10 -c --secondary=no -t 16 $B2 $W/rows.fa 2>/dev/null | awk '$1~/bet/&&$10>2000{print "B2",$1,$5,$3,$4,$6,$8,$9,$10,$11}'
echo "== F N runs near loci (deposited, +-300 kb)"
for r in "$T1 OZ418350.1:110700000-111350000" "$T2 OZ418090.1:109970000-110670000"; do set -- $r
 $S faidx $1 $2 | grep -v '>' | tr -d '\n' | python3 -c "
import sys,re
s=sys.stdin.read().upper(); st=int('$2'.split(':')[1].split('-')[0])
print('$2',len(s),[(st+m.start(),len(m.group())) for m in re.finditer('N+',s)])"
done
echo "== E reads across breakpoints"
cd gene_views/TMEM114
python3 - <<'PY'
import subprocess,re
S="/home/toki/miniforge3/envs/sinederella/bin/samtools"
out=subprocess.run([S,"view","reads.bam"],capture_output=True,text=True).stdout.splitlines()
A=[]
for l in out:
    f=l.split("\t"); ref=f[2]; p=int(f[3]); fl=int(f[1]); span=sum(int(n) for n,o in re.findall(r"(\d+)([MDN=X])",f[5]))
    cl=[int(n) for n,o in re.findall(r"(\d+)([SH])",f[5])]
    A.append((ref,p,p+span,fl,max(cl) if cl else 0,f[0]))
def sp(ref,pos,m=1000): return sum(1 for a in A if a[0].startswith(ref) and a[1]<=pos-m and a[2]>=pos+m and not a[3]&0x800)
print("hap1 junction 45709 spanning",sp("tri_hap1",45709))
print("hap2 bpL 42831",sp("tri_hap2",42831),"bpR 102848",sp("tri_hap2",102848))
for pos in range(30000,125000,10000):
    print("hap2 pos",pos,sum(1 for a in A if a[0].startswith("tri_hap2") and a[1]<=pos<a[2]),end="; ")
print()
for pos in range(30000,65000,5000):
    print("hap1 pos",pos,sum(1 for a in A if a[0].startswith("tri_hap1") and a[1]<=pos<a[2]),end="; ")
print()
PY
