#!/bin/bash --login
#SBATCH -J kallisto
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 1-0
#SBATCH -o kallisto_%j.out

# Uso (desde la carpeta que contiene trimmed/ y transcriptoma_completo/):
#   sbatch kallisto_tiempo.sh
# 1) Indexa el transcriptoma  2) Cuantifica cada muestra, una tras otra.
# Version con cronometro para medir cuanto tarda cada paso.

module load spack
module load miniconda3
source activate trinity_taller

TRANSCRIPTOME="/scratch/taller/datos_pez/transcriptoma_completo/GKNI01.1.fsa_nt.gz"
TRIMMED="trimmed"
OUT_DIR="kallisto_conteos"

mkdir -p ${OUT_DIR}

INICIO=$(date +%s)
echo "Inicio: $(date)"

# 1) Indice del transcriptoma (se hace una sola vez)
T0=$(date +%s)
kallisto index -i ${OUT_DIR}/transcriptoma.idx ${TRANSCRIPTOME}
echo "Indice terminado en $(( ($(date +%s) - T0) / 60 )) min"

# 2) Cuantificacion por muestra
for R1 in ${TRIMMED}/*_R1.trimmed.fastq.gz; do
  S=$(basename "$R1" _R1.trimmed.fastq.gz)
  echo "Procesando ${S}"
  T0=$(date +%s)

  kallisto quant -i ${OUT_DIR}/transcriptoma.idx -o ${OUT_DIR}/${S} -t 4 \
    ${TRIMMED}/${S}_R1.trimmed.fastq.gz ${TRIMMED}/${S}_R2.trimmed.fastq.gz

  echo "${S} terminado en $(( ($(date +%s) - T0) / 60 )) min"
done

TOTAL=$(( $(date +%s) - INICIO ))
echo "Fin: $(date)"
echo "Tiempo total: $(( TOTAL / 3600 )) h $(( (TOTAL % 3600) / 60 )) min"
