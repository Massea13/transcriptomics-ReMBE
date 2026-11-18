#!/bin/bash
#SBATCH -J QC_bash              # Nombre del job
#SBATCH -p short                # Partición
#SBATCH -N 1                    # Nodos
#SBATCH -c 4                    # Núcleos
#SBATCH --mem=8G                # Memoria
#SBATCH -t 24:00:00             # Tiempo máximo
#SBATCH -o QC_bash-%j.out
#SBATCH -e QC_bash-%j.err

set -euo pipefail

#Workflow para correr fastqc, trimmomatic y multiqc usando bash
#----------------------------------------
# 1. Cargar módulos del cluster (solo Java, no Nextflow)
#----------------------------------------
module load spack/1.0.2
module load jdk/21.0.2-gcc-12.2.0-linux-almalinux9-x86_64-tsyhzjs

#----------------------------------------
# 2. Activar Conda y el ambiente rembe_qc
#----------------------------------------
CONDA_DIR="$HOME/miniconda3"
if [ ! -d "$CONDA_DIR" ]; then
    echo "[ERROR] Miniconda no encontrada en $CONDA_DIR" >&2
    exit 1
fi

eval "$("$CONDA_DIR/bin/conda" shell.bash hook)"

if ! conda env list | grep -qE "^rembe_qc\s"; then
    echo "[ERROR] El ambiente 'rembe_qc' no existe. Créalo primero con 01_fastqc.sh" >&2
    exit 1
fi

conda activate rembe_qc

echo "[INFO] Entorno activo: $(which fastqc) / $(which trimmomatic) / $(which multiqc)"
fastqc --version
trimmomatic -version
multiqc --version

#----------------------------------------
# 3. Rutas del proyecto
#----------------------------------------
project="quality-assessment"
home="/home/idalia/transcriptomics-ReMBE/dia_1"
data="/scratch/taller/datos_oyamel"
scripts="/home/idalia/transcriptomics-ReMBE/dia_1/scripts"

ADAPTERS="${CONDA_PREFIX}/share/trimmomatic-0.41-0/adapters/TruSeq3-PE-2.fa"

# Verificar que los adaptadores existen
if [ ! -f "$ADAPTERS" ]; then
    echo "[ERROR] No se encontró el archivo de adaptadores: $ADAPTERS" >&2
    echo "[INFO] Búscalo con: find ${CONDA_PREFIX} -name 'TruSeq3-PE-2.fa'" >&2
    exit 1
fi

# Crear estructura de salida
RESULTS="${home}/${project}/results_bash"
mkdir -p "${RESULTS}"/{fastqc_raw,trimmed,fastqc_trimmed,multiqc,logs}
LOG_DIR="${RESULTS}/logs"

#----------------------------------------
# 4. Detectar pares R1/R2 con regex
#----------------------------------------
shopt -s nullglob
FASTQS=( "${data}/RAW/"*.fastq.gz )

if [ ${#FASTQS[@]} -eq 0 ]; then
    echo "[ERROR] No se encontraron archivos FASTQ en ${data}/RAW/" >&2
    exit 1
fi

declare -A R1 R2 SAMPLES

for f in "${FASTQS[@]}"; do
    base=$(basename "$f")
    if [[ "$base" =~ ^(.+)_R1_001\.fastq\.gz$ ]]; then
        sample="${BASH_REMATCH[1]}"
        R1["$sample"]="$f"
        SAMPLES["$sample"]=1
    elif [[ "$base" =~ ^(.+)_R2_001\.fastq\.gz$ ]]; then
        sample="${BASH_REMATCH[1]}"
        R2["$sample"]="$f"
        SAMPLES["$sample"]=1
    else
        echo "[WARN] Archivo con nombre inesperado, se ignora: $base" >&2
    fi
done

# Verificar que cada muestra tenga R1 y R2
VALID_SAMPLES=()
for s in "${!SAMPLES[@]}"; do
    if [[ -n "${R1[$s]:-}" && -n "${R2[$s]:-}" ]]; then
        VALID_SAMPLES+=("$s")
    else
        echo "[WARN] Muestra $s incompleta (R1=${R1[$s]:-NO}, R2=${R2[$s]:-NO}); se omite" >&2
    fi
done

if [ ${#VALID_SAMPLES[@]} -eq 0 ]; then
    echo "[ERROR] No hay muestras con R1 y R2 completos" >&2
    exit 1
fi

# Ordenar muestras para reproducibilidad
IFS=$'\n' VALID_SAMPLES=($(sort <<<"${VALID_SAMPLES[*]}")); unset IFS

echo "[INFO] Muestras válidas (${#VALID_SAMPLES[@]}):"
printf '  - %s\n' "${VALID_SAMPLES[@]}"

#----------------------------------------
# 5. Paso 1: FastQC sobre raw
#----------------------------------------
echo "[STEP 1] FastQC raw ..."
for s in "${VALID_SAMPLES[@]}"; do
    echo "  -> $s"
    fastqc --threads 4 --outdir "${RESULTS}/fastqc_raw" \
        "${R1[$s]}" "${R2[$s]}" \
        > "${LOG_DIR}/${s}.fastqc_raw.log" 2>&1
done

#----------------------------------------
# 6. Paso 2: Trimmomatic
#----------------------------------------
echo "[STEP 2] Trimmomatic ..."
for s in "${VALID_SAMPLES[@]}"; do
    echo "  -> $s"
    trimmomatic PE -threads 4 -phred33 \
        "${R1[$s]}" "${R2[$s]}" \
        "${RESULTS}/trimmed/${s}_R1_paired.fastq.gz" \
        "${RESULTS}/trimmed/${s}_R1_unpaired.fastq.gz" \
        "${RESULTS}/trimmed/${s}_R2_paired.fastq.gz" \
        "${RESULTS}/trimmed/${s}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
        SLIDINGWINDOW:4:20 \
        MINLEN:36 \
        > "${LOG_DIR}/${s}.trimmomatic.log" 2>&1
done

#----------------------------------------
# 7. Paso 3: FastQC sobre trimmed (solo paired)
#----------------------------------------
echo "[STEP 3] FastQC trimmed ..."
for s in "${VALID_SAMPLES[@]}"; do
    echo "  -> $s"
    fastqc --threads 4 --outdir "${RESULTS}/fastqc_trimmed" \
        "${RESULTS}/trimmed/${s}_R1_paired.fastq.gz" \
        "${RESULTS}/trimmed/${s}_R2_paired.fastq.gz" \
        > "${LOG_DIR}/${s}.fastqc_trimmed.log" 2>&1
done

#----------------------------------------
# 8. Paso 4: MultiQC (raw + trimmed)
#----------------------------------------
echo "[STEP 4] MultiQC ..."
multiqc \
    "${RESULTS}/fastqc_raw" \
    "${RESULTS}/fastqc_trimmed" \
    -o "${RESULTS}/multiqc" \
    -n multiqc_report.html \
    > "${LOG_DIR}/multiqc.log" 2>&1

echo "[INFO] Quality assessment completado."
echo "[INFO] Resultados en: ${RESULTS}"
echo "[INFO] Informe MultiQC: ${RESULTS}/multiqc/multiqc_report.html"
