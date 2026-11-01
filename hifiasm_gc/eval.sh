#!/bin/bash
# Evaluate GC-corrected hifiasm haplotypes (asm_gc) vs stock asm25.
set -euo pipefail
G=/home/toki/gof
O=$G/hifiasm_gc/eval
A=$G/asm25
N=$G/hifiasm_gc/asm_gc
C=$G/code
mkdir -p $O
cd $O || exit 1
export PATH=/home/toki/miniforge3/envs/asm25/bin:$PATH:/home/toki/miniforge3/envs/hifiasm/bin
which minimap2 compleasm merqury.sh
export MERQURY=$(ls -d /home/toki/miniforge3/envs/asm25/share/merqury*)
export TMPDIR=/home/toki/tmp
T=${1:-24}
step() { echo "#### $* -- $(date '+%F %T')"; }
step "sizes"
for f in $A/tri.hic.hap1.p_ctg.fa $A/tri.hic.hap2.p_ctg.fa $N/tri.hic.hap1.p_ctg.fa $N/tri.hic.hap2.p_ctg.fa; do
  awk -v f=$f '/^>/{n++; next}{L+=length($0); a[n]+=length($0)} END{print f, n, L}' $f
done | tee sizes.txt
step "unitig rule test"
[ -s utg_rule_gc.tsv ] || python3 $C/hifi/utg_rule_test.py $N/tri.hic.p_utg.gfa $N/tri.hic.p_utg.noseq.gfa $N/tri.hic.hap1.p_ctg.noseq.gfa $N/tri.hic.hap2.p_ctg.noseq.gfa > utg_rule_gc.tsv
[ -s utg_rule_asm25.tsv ] || python3 $C/hifi/utg_rule_test.py $A/tri.hic.p_utg.gfa $A/tri.hic.p_utg.noseq.gfa $A/tri.hic.hap1.p_ctg.noseq.gfa $A/tri.hic.hap2.p_ctg.noseq.gfa > utg_rule_asm25.tsv
cmp $N/tri.hic.p_utg.noseq.gfa $A/tri.hic.p_utg.noseq.gfa && echo "p_utg identical to asm25"
step "10 regions + controls: windows vs both haplotypes"
for h in hap1 hap2; do
  minimap2 -x asm5 -t $T $A/tri.hic.$h.p_ctg.fa $G/ends/sv/windows.fa > win_asm25_$h.paf 2>> mm2.log
  minimap2 -x asm5 -t $T $N/tri.hic.$h.p_ctg.fa $G/ends/sv/windows.fa > win_gc_$h.paf 2>> mm2.log
done
python3 $G/hifiasm_gc/winpaf.py asm25_h1=win_asm25_hap1.paf asm25_h2=win_asm25_hap2.paf gc_h1=win_gc_hap1.paf gc_h2=win_gc_hap2.paf > windows.tsv
cat windows.tsv
step "gapfill of deposited gaps with asm25 and GC haplotypes"
for d in dep1 dep2; do
  if [ $d = dep1 ]; then R=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna; else R=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna; fi
  mkdir -p gf; cd $O/gf || exit 1
  python3 $C/wgs/gapfill.py $R $d asm_h1=$A/tri.hic.hap1.p_ctg.fa asm_h2=$A/tri.hic.hap2.p_ctg.fa gc_h1=$N/tri.hic.hap1.p_ctg.fa gc_h2=$N/tri.hic.hap2.p_ctg.fa --threads $T > $d.log 2>&1
  cd $O || exit 1
  python3 $G/hifiasm_gc/gapcmp.py gf/$d/gaps.tsv $G/gf_all/flankcn/gap_causes.tsv $G/repair/repairs.tsv $d > gapcmp_$d.txt
  head -20 gapcmp_$d.txt
done
step "compleasm (contigs)"
for x in gc:$N; do
  n=${x%%:*}; p=${x#*:}
  for h in hap1 hap2; do
    [ -s cpl_${n}_$h/summary.txt ] || compleasm run -a $p/tri.hic.$h.p_ctg.fa -o cpl_${n}_$h -l mammalia --odb odb10 -L $G/compleasm_lib -t $T > cpl_${n}_$h.log 2>&1
    echo "$n $h $(grep -E '^[SDFM]:' cpl_${n}_$h/summary.txt | tr '\n' ' ')"
  done
done | tee compleasm.txt
step "merqury (both haplotypes together)"
for x in gc:$N; do
  n=${x%%:*}; p=${x#*:}
  mkdir -p mq_$n; cd $O/mq_$n || exit 1
  ln -sf $p/tri.hic.hap1.p_ctg.fa ${n}_h1.fa; ln -sf $p/tri.hic.hap2.p_ctg.fa ${n}_h2.fa
  [ -s mq.qv ] || OMP_NUM_THREADS=$T merqury.sh $G/sqc/hifi.meryl ${n}_h1.fa ${n}_h2.fa mq > merqury.log 2>&1 || echo "merqury $n exit $?"
  cat mq.completeness.stats mq.qv || true
  cd $O || exit 1
done
step done
