#!/bin/bash --login
#SBATCH -J trimmomatic
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 1-0
#SBATCH -o trimmomatic_%j.out

module load spack
module load trimmomatic

RAW_READS="/scratch/taller/datos_pez/datos_crudos_pez"
OUT_DIR="trimmed"
ADAPTERS="/opt/ohpc/pub/apps/spack/local/linux-almalinux9-haswell/gcc-12.2.0/trimmomatic-0.39-cfbxtlbjqkajus4ezyeocgbkqxm>

mkdir -p ${OUT_DIR}

for R1 in ${RAW_READS}/*_1.fastq.gz; do
  S=$(basename "$R1" _1.fastq.gz)

  trimmomatic PE -phred33 -threads 4 \
    ${RAW_READS}/${S}_1.fastq.gz ${RAW_READS}/${S}_2.fastq.gz \
    ${OUT_DIR}/${S}_R1.trimmed.fastq.gz ${OUT_DIR}/${S}_R1.unpaired.fastq.gz \
    ${OUT_DIR}/${S}_R2.trimmed.fastq.gz ${OUT_DIR}/${S}_R2.unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50
done
