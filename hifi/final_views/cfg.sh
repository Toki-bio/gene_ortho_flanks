# paths for the final-assembly views
G=/home/toki/gof; W=$G/final_views; F=$G/final; LR=$G/lr_all; GL=$G/geneloss3
T1=/home/toki/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
T2=/home/toki/Sicista2026/demog/genomes/tri_hap2.fna
BET=/home/toki/Sicista2026/assembly/primary_eval/hifiasm_primary.p_ctg.fa
READS=$(ls $G/reads/*.fastq.gz | head -1)
SCF=/home/toki/sineclose_bash/fast3/scfish2
export LC_ALL=C
export PATH=/home/toki/miniforge3/envs/sinederella/bin:/home/toki/miniforge3/envs/hifiasm/bin:/home/toki/miniforge3/envs/jf/bin:/home/toki/miniforge3/envs/rapidgz/bin:$PATH
export SCFISH_DEC="rapidgzip -d -c -P 16"
FL=25000
