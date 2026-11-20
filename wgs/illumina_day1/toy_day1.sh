#!/bin/bash
# Toy test of illumina_day1.sh on one betulina contig with simulated diploid Illumina reads.
# Truth: hap A = contig, hap B = contig + the long-read het SNVs; 30x 2x150 wgsim reads (15x per haplotype, 0.2% errors, no indels).
# The toy ASSEMBLY gets 20 planted 1-bp deletions in exon homopolymers (5-8 bp) in the last third of the contig.
# Expected: same-animal 'both' high; QV before < after; the 20 planted tracts called 'assembly wrong' (+1) and fixed by pypolca.
set -euo pipefail
D=$(dirname $(readlink -f $0))
C=${1:-ptg000170l}; W=/home/toki/Sicista2026/illumina_prep/toy
mkdir -p $W; cd $W || exit 1
export PATH=/home/toki/miniforge3/envs/polish/bin:$PATH TMPDIR=/home/toki/tmp
PY=/home/toki/Sicista2026/demog/envs/smcpp/bin/python
P=/home/toki/Sicista2026/assembly/primary_eval/hifiasm_primary.p_ctg.fa
samtools faidx $P $C > hapA.fa
zcat /home/toki/Sicista2026/demog/bet_single/smcpp/het.vcf.gz | awk -v c=$C '$1==c{print $1"\t"$2"\t"$4"\t"$5}' > hets.tsv
gzip -c hets.tsv > sites.tsv.gz
awk -v c=$C '$1==c' /home/toki/Sicista2026/hiller_share/TOGA2/query_annotation.bed > exons.bed
awk -v c=$C 'NR==1 || $1==c' /home/toki/Sicista2026/demog/hp/bet.classes.tsv > lr_classes.tsv
$PY - <<EOF
import random
random.seed(5)
seq=''.join(l.strip() for l in open('hapA.fa') if l[0]!='>'); L=len(seq)
b=list(seq)
for l in open('hets.tsv'):
    c,p,r,a=l.split(); p=int(p)-1
    if b[p].upper()==r: b[p]=a
open('hapB.fa','w').write('>$C\n'+''.join(b)+'\n')
# exon homopolymers 5-8 bp in the last third, from the long-read tract table (keys are positions in the unmodified contig)
cand=[l.split('\t') for l in open('lr_classes.tsv')][1:]
cand=[(int(f[1]),int(f[3])) for f in cand if int(f[1])>2*L//3 and 5<=int(f[3])<=8]
pick=sorted(random.sample(cand,20))
s=list(seq)
for p,n in sorted(pick,reverse=True): del s[p]          # shorten each tract by one base (deleting its first base)
open('toy_asm.fa','w').write('>$C\n'+''.join(s)+'\n')
with open('planted.tsv','w') as o:
    sh=0
    for p,n in pick: o.write(f'{p}\t{p-sh}\t{n}\n'); sh+=1        # original pos0, pos0 in toy assembly, tract length
print('planted',len(pick),'contig length',L, 'hets', sum(1 for _ in open('hets.tsv')))
EOF
# exons BED in toy coordinates: shift blocks after each deletion (planted deletions are 1 bp each)
$PY - <<EOF
import bisect
pl=sorted(int(l.split()[0]) for l in open('planted.tsv'))
def sh(x): return x-bisect.bisect_left(pl,x)
o=open('exons_toy.bed','w')
for l in open('exons.bed'):
    f=l.rstrip('\n').split('\t'); cs,ce=int(f[1]),int(f[2])
    st=[int(x) for x in f[11].strip(',').split(',')]; sz=[int(x) for x in f[10].strip(',').split(',')]
    a=[sh(cs+s) for s in st]; e=[sh(cs+s+z) for s,z in zip(st,sz)]
    f[1]=str(a[0]); f[2]=str(e[-1]); f[6]=str(sh(int(f[6]))); f[7]=str(sh(int(f[7])))
    f[10]=','.join(str(y-x) for x,y in zip(a,e))+','; f[11]=','.join(str(x-a[0]) for x in a)+','
    o.write('\t'.join(f)+'\n')
EOF
# shift the het sites and the long-read tract table into toy-assembly coordinates (1-bp deletions at planted.tsv col 1)
$PY - <<'EOF2'
import bisect, gzip
pl=sorted(int(l.split()[0]) for l in open('planted.tsv'))
def sh(x): return x-bisect.bisect_left(pl,x)
with gzip.open('sites_toy.tsv.gz','wt') as o:
    for l in open('hets.tsv'):
        c,p,r,a=l.split(); o.write(f'{c}\t{sh(int(p)-1)+1}\t{r}\t{a}\n')
with open('lr_classes_toy.tsv','w') as o:
    for i,l in enumerate(open('lr_classes.tsv')):
        f=l.rstrip('\n').split('\t')
        if i: f[1]=str(sh(int(f[1])))
        o.write('\t'.join(f)+'\n')
EOF2
L=$(awk '!/^>/{n+=length($0)} END{print n}' hapA.fa); N=$((L*15/300))
wgsim -N $N -1 150 -2 150 -e 0.002 -r 0 -R 0 -X 0 -S 11 hapA.fa a_1.fq a_2.fq > /dev/null
wgsim -N $N -1 150 -2 150 -e 0.002 -r 0 -R 0 -X 0 -S 12 hapB.fa b_1.fq b_2.fq > /dev/null
cat a_1.fq b_1.fq | gzip > r1.fq.gz; cat a_2.fq b_2.fq | gzip > r2.fq.gz; rm -f a_?.fq b_?.fq
samtools faidx toy_asm.fa; bwa-mem2 index toy_asm.fa > /dev/null 2>&1
bash $D/illumina_day1.sh r1.fq.gz r2.fq.gz out 16 toy_asm.fa sites_toy.tsv.gz lr_classes_toy.tsv exons_toy.bed > day1.log 2>&1 || { tail -30 day1.log; exit 1; }
tail -40 day1.log
echo "== planted tracts in the Illumina table (toy coordinates)"
awk 'NR==FNR{p[$2]=1; next} FNR>1 && ($2 in p){print $2, $4, $5, $6, $9, $10}' planted.tsv out/ill_tracts.classes.tsv | sort -n | head -25
echo "== planted sites in the polished assembly: tract restored?"
$PY - <<EOF
import pysam
pa=pysam.FastaFile('out/polished.fa'); ta=pysam.FastaFile('toy_asm.fa'); h=pysam.FastaFile('hapA.fa')
c=pa.references[0]; ok=0
for l in open('planted.tsv'):
    o,t,n=map(int,l.split())
    truth=h.fetch(c,o-10,o+int(n)+10)
    ok += truth in pa.fetch(c,max(0,t-2000),t+2000)
print('planted tracts restored by pypolca:',ok,'of 20')
EOF
