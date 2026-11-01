#!/bin/bash
# LJA 0.2 (bioconda, env /home/toki/miniforge3/envs/lja) toys on the round-1 fished reads. Result 2026-11-01: all three
# runs (--diploid fastq.gz, haploid fastq.gz, --diploid FASTA) end with "Child process crashed" in "Uncompressing
# homopolymers in contigs" after ~3-5 min; only the homopolymer-compressed mdbg/assembly.hpc.fasta is written, which
# cannot serve as fill sequence. LJA was therefore not used; HiCanu (asm2_hicanu.sh) is the independent assembler.
set -euo pipefail
cd /home/toki/gof/asm2/toy || exit 1
export TMPDIR=/home/toki/tmp
L=/home/toki/miniforge3/envs/lja/bin/lja
[ -s fished.fa ] || zcat /home/toki/gof/repair/fished.fq.gz | awk 'NR%4==1{print ">" substr($1,2)} NR%4==2{print}' > fished.fa
/usr/bin/time -v $L --diploid -t 32 -o lja --reads /home/toki/gof/repair/fished.fq.gz > lja.out 2> lja.time || true
/usr/bin/time -v $L -t 16 -o hap --reads /home/toki/gof/repair/fished.fq.gz > hap.out 2> hap.time || true
/usr/bin/time -v $L --diploid -t 16 -o dipfa --reads fished.fa > dipfa.out 2> dipfa.time || true
for x in lja hap dipfa; do echo "== $x"; tail -n 2 $x.out; done
