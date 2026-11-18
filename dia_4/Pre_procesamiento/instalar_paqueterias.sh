#!/bin/bash --login
#SBATCH -J instalar_R
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 0-2
#SBATCH -o instalar_R_%j.out

# Uso (una sola vez, desde la carpeta que contiene 00_instalar_paquetes.R):
#   sbatch instalar_paquetes.sh
# Al terminar, revisa el final de instalar_R_<numero>.out: debe decir "Todo listo."

module swap gnu12 gnu14/14.2.0
module load R/4.4.1

# Carpeta personal para los paquetes de R (la del sistema no se puede modificar).
# Tiene que existir; si no, R la ignora.
export R_LIBS_USER=~/R/4.4
mkdir -p ${R_LIBS_USER}

# Compila usando los 4 nucleos (mas rapido)
export MAKEFLAGS="-j4"

Rscript 00_instalar_paqueterias.R