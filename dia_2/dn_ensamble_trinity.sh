#!/bin/bash --login
#SBATCH -p multicore
#SBATCH -J trinity_denovo
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH -t 2-0
#SBATCH -o trinity_denovo_%j.out
#SBATCH -e trinity_denovo_%j.err

# Ensamble de novo de un transcriptoma con Trinity a partir de lecturas
# paired-end ya limpias (sin adaptadores y con recorte por calidad).
# Ajusta rutas y parametros segun tus datos y tu cluster.

# -------------------------
# Módulos
# -------------------------
module load apps/singularity/trinity/2.15.1

# -------------------------
# paths y parámetros
# -------------------------
READS_DIR="trimmed_reads"
R1_SUFFIX="_R1.trimmed.fastq.gz"
R2_SUFFIX="_R2.trimmed.fastq.gz"

# El nombre debe incluir "trinity"; de lo contrario, Trinity no se ejecuta
OUTPUT_DIR="trinity_denovo"

THREADS="${SLURM_CPUS_PER_TASK}"
MAX_MEMORY="90G"          # Menor que --mem para dejar un margen
MIN_CONTIG_LENGTH=200     # Longitud minima de los contigs reportados
MAX_READ_COV=200          # Cobertura maxima tras la normalizacion in silico

# Orientacion de la libreria:
#   dUTP / Illumina TruSeq Stranded: STRAND_OPT="--SS_lib_type RF"
#   Libreria no direccional:         STRAND_OPT=""
STRAND_OPT="--SS_lib_type RF"

# -------------------------
# Listas de lecturas
# -------------------------
# Trinity recibe las lecturas de todas las muestras en dos listas separadas
# por comas. Ambas deben conservar el mismo orden de muestras.
shopt -s nullglob
R1_FILES=( "${READS_DIR}"/*"${R1_SUFFIX}" )
shopt -u nullglob

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
  echo "ERROR: no hay archivos *${R1_SUFFIX} en ${READS_DIR}" >&2
  exit 1
fi

R2_FILES=()
for R1 in "${R1_FILES[@]}"; do
  R2="${R1%"${R1_SUFFIX}"}${R2_SUFFIX}"
  if [[ ! -f "${R2}" ]]; then
    echo "ERROR: falta el par de ${R1}: ${R2}" >&2
    exit 1
  fi
  R2_FILES+=( "${R2}" )
done

LEFT=$(IFS=,; echo "${R1_FILES[*]}")
RIGHT=$(IFS=,; echo "${R2_FILES[*]}")

echo "Muestras: ${#R1_FILES[@]}"
echo "Left:  ${LEFT}"
echo "Right: ${RIGHT}"

# -------------------------
# Ensamble de novo
# -------------------------
trinity Trinity \
  --seqType fq \
  --left "${LEFT}" \
  --right "${RIGHT}" \
  ${STRAND_OPT} \
  --min_contig_length "${MIN_CONTIG_LENGTH}" \
  --normalize_max_read_cov "${MAX_READ_COV}" \
  --max_memory "${MAX_MEMORY}" \
  --CPU "${THREADS}" \
  --output "${OUTPUT_DIR}"