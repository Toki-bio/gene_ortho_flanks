#!/bin/bash
# Genome-wide runs: depth-vs-GC inputs for betulina (ONT, existing BAM) and trizona (HiFi, all reads -> deposited hap1),
# trizona hap1-hap2 SNP bins (two-copy proxy), betulina primary vs hap1/hap2 (noscaf), betulina contig-end inputs,
# and flank mapping of trizona gaps / random trizona sites to the three betulina assemblies.
set -euo pipefail
cd /home/toki/gof/betulina || exit 1
E=/home/toki/miniforge3/envs
export PATH=$E/sinederella/bin:$PATH
MM=$E/hifiasm/bin/minimap2
MOS=$E/lrtools/bin/mosdepth
C=/home/toki/gof/betulina/code
HC=/home/toki/gof/code/hifi
P=/home/toki/Sicista2026/assembly/primary_eval/hifiasm_primary.p_ctg.fa
H1=/home/toki/Sicista2026/assembly/hifiasm_out/sicista_noscaf.p_ctg.fa
H2GFA=/home/toki/Sicista2026/assembly/hifiasm_out/sicista_noscaf.bp.hap2.p_ctg.gfa
BAM=/home/toki/Sicista2026/qc/sbet_ont.bam
VCF=/home/toki/Sicista2026/qc/calls_full/clair_auto/merge_output.vcf.gz
T1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
T2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
HIFI=/home/toki/gof/reads/ERR16751478.fastq.gz
TGAPS=/home/toki/gof/gf_all/dep1/gaps.tsv
TFL=/home/toki/gof/gf_all/flankcn/flanks.fa
mkdir -p gw && cd gw || exit 1
date
# A: trizona HiFi depth, all reads (16 threads)
( zcat $HIFI | $MM -x map-hifi -t 16 --secondary=no $T1 - 2> tri_hifi.mm.log \
   | python3 $C/paf_depth.py $T1.fai 10000 20 > tri_hifi.depth10k.tmp && mv tri_hifi.depth10k.tmp tri_hifi.depth10k.bed && echo A_ok ) &
PA=$!
# B: trizona hap1 vs hap2 SNPs per hap1 bin (4 threads)
( python3 $C/chunk_fa.py $T2 1000000 | $MM -cx asm5 --cs -t 6 $T1 - 2> tri_hap.mm.log | python3 $C/hap_snp_bins.py 10000 > tri.hapsnp.tsv && echo B_ok ) &
PB=$!
# C: betulina ONT depth (MAPQ>=20, 10 kb) + Clair3 het bins (4 threads)
( $MOS -n -t 4 -Q 20 --by 10000 bet $BAM && zcat bet.regions.bed.gz > bet.depth10k.bed \
   && zcat $VCF | python3 $C/het_bins.py 10000 > bet.het.tsv && echo C_ok ) &
PC=$!
python3 $C/bins_gc.py $P 10000 > bet.gc.bed
python3 $C/bins_gc.py $T1 10000 > tri.gc.bed
# D: betulina hap2 (noscaf) fasta; primary vs hap1 / hap2 (2 x 4 threads)
[ -s bet_hap2.fa ] || awk '$1=="S"{print ">"$2; print $3}' $H2GFA > bet_hap2.fa
samtools faidx bet_hap2.fa
( python3 $C/chunk_fa.py $P 1000000 | $MM -x asm5 -t 4 --secondary=no $H1 - 2> p_h1.mm.log > p_vs_h1.paf && echo D1_ok ) &
PD1=$!
( python3 $C/chunk_fa.py $P 1000000 | $MM -x asm5 -t 4 --secondary=no bet_hap2.fa - 2> p_h2.mm.log > p_vs_h2.paf && echo D2_ok ) &
PD2=$!
# E: contig-end inputs (primary) - cheap
python3 $C/ends_prep.py $P ends 30000 20000
samtools bedcov -Q 20 ends.bins.bed $BAM > ends.bedcov
samtools view -q 20 -M -L ends.reads3k.bed $BAM | python3 $C/ends_reads.py $P.fai > ends.reads.tsv
trf ends.ends.fa 2 7 7 80 10 50 500 -h -ngs > ends.trf.dat || true
# random trizona control sites + flanks
python3 $C/random_sites.py $T1.fai $TGAPS 1000 7 > rnd.tsv
python3 $HC/gap_flank_cn.py prep $T1 rnd.tsv rnd_flanks.fa
cat $TFL rnd_flanks.fa ends.ends.fa > q_all.fa
wait $PB $PC $PD1 $PD2 || echo SOME_BG_FAILED
date; echo BCD_done
# F: flanks (trizona gaps, random sites, betulina ends) -> betulina primary / hap1 / hap2; betulina ends -> trizona hap1
$MM -x asm20 -t 12 $P q_all.fa 2> qa_p.mm.log > qa_vs_p.paf
$MM -x asm20 -t 12 $H1 q_all.fa 2> qa_h1.mm.log > qa_vs_h1.paf
$MM -x asm20 -t 12 bet_hap2.fa q_all.fa 2> qa_h2.mm.log > qa_vs_h2.paf
$MM -x asm20 -t 12 $T1 ends.ends.fa 2> ends_t1.mm.log > ends_vs_t1.paf
date; echo F_done
wait $PA || echo A_FAILED
date; echo ALL_DONE
