#!/bin/bash --login
#SBATCH -J kallisto
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 1-0
#SBATCH -o kallisto_%j.out

# Cuantificación de transcritos con kallisto
# Uso (desde la carpeta que contiene trimmed/):
#   sbatch kallisto.sh
#
# 1) Indexa el transcriptoma
# 2) Cuantifica cada muestra pareada (R1/R2)

module load spack
module load miniconda3
source activate trinity_taller

# --- Parámetros ---
TRANSCRIPTOME="/scratch/taller/datos_pez/transcriptoma_completo/GKNI01.1.fsa_nt.gz"
TRIMMED="/scratch/taller/datos_pez/trimmed"
OUT_DIR="kallisto_conteos"
INDEX="${OUT_DIR}/transcriptoma.idx"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "${OUT_DIR}"

# --- 1) Índice del transcriptoma ---
kallisto index -i "${INDEX}" "${TRANSCRIPTOME}"

# --- 2) Cuantificación por muestra ---
for R1 in "${TRIMMED}"/*_R1.trimmed.fastq.gz; do
  S=$(basename "${R1}" _R1.trimmed.fastq.gz)
  R2="${TRIMMED}/${S}_R2.trimmed.fastq.gz"

  echo "Procesando ${S}"
  kallisto quant -i "${INDEX}" -o "${OUT_DIR}/${S}" -t "${THREADS}" "${R1}" "${R2}"
done

echo "Listo. Resultados en ${OUT_DIR}/"
