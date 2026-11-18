#!/bin/bash --login
#SBATCH -J R_analisis
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 0-1
#SBATCH -o R_analisis_%j.out

module swap gnu12 gnu14/14.2.0
module load R/4.4.1

export R_LIBS_USER=~/R/4.4
mkdir -p ${R_LIBS_USER}

Rscript 00_exportar_counts.R