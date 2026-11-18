#!/bin/bash --login
#SBATCH -J R_analisis
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 0-1
#SBATCH -o R_analisis_%j.out

# Uso (desde la carpeta que contiene kallisto_counts/ y samples.tsv):
#   sbatch R_analisis.sh
# Los paquetes se instalan antes, una sola vez, con 00_instalar_paquetes.R

module swap gnu12 gnu14/14.2.0
module load R/4.4.1
export R_LIBS_USER=~/R/4.4

Rscript ../00_instalar_paquetes.R
Rscript ../00_exportar_counts.R