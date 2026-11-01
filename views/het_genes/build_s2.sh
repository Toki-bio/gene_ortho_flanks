set -euo pipefail
cd /home/toki/gof/views_het || exit 1
E=/home/toki/miniforge3/envs
SAM=$E/sinederella/bin/samtools; MM=$E/hifiasm/bin/minimap2; MAFFT=$(ls $E/sinederella/bin/mafft $E/*/bin/mafft 2>/dev/null | head -1)
T1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
B1=/home/toki/Sicista2026/assembly/hifiasm_out/sicista_hifiasm.p_ctg.fa
B2=/home/toki/gof/betulina/gw/bet_hap2.fa
[ -f $B1.fai ] || $SAM faidx $B1
rc(){ python3 -c "
import sys
s=[l.strip() for l in sys.stdin if not l.startswith('>')]
s=''.join(s).upper(); t=s[::-1].translate(str.maketrans('ACGTN','TGCAN'))
print('>'+sys.argv[1]); [print(t[i:i+80]) for i in range(0,len(t),80)]" "$1"; }
fa(){ python3 -c "
import sys
s=''.join(l.strip() for l in sys.stdin if not l.startswith('>')).upper()
print('>'+sys.argv[1]); [print(s[i:i+80]) for i in range(0,len(s),80)]" "$1"; }
mk(){ # gene t2file hap1wide b1reg b1strand b2reg b2strand
 g=$1
 fa "${g}_tri_hap2" < $2 > $g.r2.fa
 $SAM faidx $T1 $3 > $g.h1wide.fa
 $MM -x asm20 -c --secondary=no $g.r2.fa $g.h1wide.fa 2>/dev/null | cut -f1-12 | awk '$10>300' > $g.h1.paf
 cat $g.h1.paf
 python3 - $g <<'PY'
import sys
g=sys.argv[1]
rows=[l.split('\t') for l in open(g+'.h1.paf')]
# strand of the majority
plus=sum(int(r[9]) for r in rows if r[4]=='+'); minus=sum(int(r[9]) for r in rows if r[4]=='-')
st='+' if plus>=minus else '-'
rows=[r for r in rows if r[4]==st]
s=min(int(r[2]) for r in rows); e=max(int(r[3]) for r in rows)
open(g+'.h1.cut','w').write(f"{st} {s} {e}\n"); print(g,'hap1 strand',st,'query-span',s,e,'wide len',rows[0][1])
PY
 read st s e < $g.h1.cut
 python3 - $g $st $s $e <<'PY'
import sys
g,st,s,e=sys.argv[1:5]; s=int(s); e=int(e)
seq=''.join(l.strip() for l in open(g+'.h1wide.fa') if not l.startswith('>')).upper()
# need target-side coords of the span: use first/last paf rows
rows=[l.split('\t') for l in open(g+'.h1.paf') if l.split('\t')[4]==st]
qs=min(int(r[2]) for r in rows); qe=max(int(r[3]) for r in rows)
seq=seq[qs:qe]
if st=='-': seq=seq[::-1].translate(str.maketrans('ACGTN','TGCAN'))
hdr=open(g+'.h1wide.fa').readline().strip()[1:]
open(g+'.r1.fa','w').write('>'+g+'_tri_hap1\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')
open(g+'.r1.info','w').write(f"{hdr} q{qs}-{qe} strand {st} len {len(seq)}\n")
PY
 cat $g.r1.info
 $SAM faidx $B1 $4 | fa "${g}_bet_hap1" > $g.r3.fa
 if [ "$5" = "-" ]; then $SAM faidx $B1 $4 | rc "${g}_bet_hap1" > $g.r3.fa; fi
 $SAM faidx $B2 $6 | fa "${g}_bet_hap2" > $g.r4.fa
 if [ "$7" = "-" ]; then $SAM faidx $B2 $6 | rc "${g}_bet_hap2" > $g.r4.fa; fi
 cat $g.r1.fa $g.r2.fa $g.r3.fa $g.r4.fa > $g.in.fa
 grep -c '>' $g.in.fa
 $MAFFT --thread 8 --auto --quiet $g.in.fa > $g.aln.raw.fa
 grep '>' $g.aln.raw.fa
}
mk TMEM114 tmem_t2.fa OZ418350.1:111007000-111051000 h1tg000312l:6327125-6423689 + h2tg001095l:630000-722000 -
mk SMCO1 smco_t2.fa OZ418345.1:165282000-165325000 h1tg000009l:16644000-16685000 + h2tg000132l:8114500-8157000 -
