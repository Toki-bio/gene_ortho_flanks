#!/bin/bash
# After main1.sh: depth-vs-GC curves, contig-end classification, primary-vs-haplotype bins, trizona sites on betulina.
set -euo pipefail
cd /home/toki/gof/betulina/gw || exit 1
E=/home/toki/miniforge3/envs
export PATH=$E/sinederella/bin:$PATH
C=/home/toki/gof/betulina/code
P=/home/toki/Sicista2026/assembly/primary_eval/hifiasm_primary.p_ctg.fa
H1=/home/toki/Sicista2026/assembly/hifiasm_out/sicista_noscaf.p_ctg.fa
BAM=/home/toki/Sicista2026/qc/sbet_ont.bam
T1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
TGAPS=/home/toki/gof/gf_all/dep1/gaps.tsv
TCAUSE=/home/toki/gof/gf_all/flankcn/gap_causes.tsv
date
# depth vs GC (betulina: exclude half-depth contigs = X; trizona: exclude chrX OZ418353.1)
python3 $C/gc_curve.py bet.depth10k.bed bet.gc.bed bet.het.tsv clair /home/toki/Sicista2026/qc/half_ge1M.txt > bet.curve.txt
echo OZ418353.1 > tri_x.txt
[ -s tri_hifi.depth10k.bed ] && python3 $C/gc_curve.py tri_hifi.depth10k.bed tri.gc.bed tri.hapsnp.tsv hapsnp tri_x.txt > tri_hifi.curve.txt
# contig ends
python3 $C/ends_classify.py ends.ends.fa ends.bedcov bet.curve.txt bet.het.tsv bet.depth10k.bed ends.trf.dat \
  ends.reads.tsv qa_vs_h1.paf qa_vs_h2.paf ends_vs_t1.paf $T1.fai $TGAPS $TCAUSE > end_causes.tsv 2> end_causes.counts
# primary vs hap1 / hap2
python3 $C/pvh.py p_vs_h1.paf p_vs_h2.paf bet.depth10k.bed bet.gc.bed bet.het.tsv bet.curve.txt > pvh.txt
# trizona gap sites + random sites -> betulina primary / hap1 / hap2
( echo site; cut -f1 $TCAUSE | tail -n +2; cut -f1 rnd.tsv | tail -n +2 ) > sites.txt
python3 $C/orth_sites.py qa_vs_p.paf $P.fai sites.txt > sites_P.tsv
python3 $C/orth_sites.py qa_vs_h1.paf $H1.fai sites.txt > sites_H1.tsv
python3 $C/orth_sites.py qa_vs_h2.paf bet_hap2.fa.fai sites.txt > sites_H2.tsv
python3 $C/site_depth.py prep sites_P.tsv $P site_bins.bed
samtools bedcov -Q 20 site_bins.bed $BAM > site_bins.bedcov
python3 $C/site_depth.py summary site_bins.bed site_bins.bedcov bet.curve.txt $TCAUSE sites_P.tsv sites_H1.tsv sites_H2.tsv > site_summary.tsv
date; echo MAIN2_DONE
