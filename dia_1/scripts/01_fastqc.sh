#!/bin/bash
#SBATCH -J QC
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 24:00:00
#SBATCH -o QC-%j.out
#SBATCH -e QC-%j.err

set -euo pipefail

#----------------------------------------
# 1. Cargar módulos del cluster (sin duplicar spack ni cargar miniconda3 del cluster)
#----------------------------------------
module load spack/1.0.2
module load jdk/21.0.2-gcc-12.2.0-linux-almalinux9-x86_64-tsyhzjs
module load nextflow

#----------------------------------------
# 2. Cargar el conda base del usuario (solo para tener el comando `conda`)
#----------------------------------------
CONDA_DIR="$HOME/miniconda3"

if [ ! -d "$CONDA_DIR" ]; then
    echo "[ERROR] No se encontró Miniconda en $CONDA_DIR." >&2
    echo "[INFO] Instálala o carga el módulo de conda del cluster." >&2
    exit 1
fi

eval "$("$CONDA_DIR/bin/conda" shell.bash hook)"

#----------------------------------------
# 3. Activar el entorno COMPARTIDO por ruta absoluta
#----------------------------------------
SHARED_ENV="/home/idalia/miniconda3/envs/rembe_qc"

if [ ! -d "$SHARED_ENV" ]; then
    echo "[ERROR] El entorno compartido no existe: $SHARED_ENV" >&2
    exit 1
fi

conda activate "$SHARED_ENV"

# Verificar que las herramientas están disponibles
fastqc --version
multiqc --version
trimmomatic -version

#----------------------------------------
# 4. Forzar Java del módulo (sobrescribe JAVA_HOME del entorno compartido)
#----------------------------------------
export JAVA_HOME=$(dirname $(dirname $(readlink -f $(which java))))
export JAVA_CMD="${JAVA_HOME}/bin/java"
export PATH="${JAVA_HOME}/bin:${PATH}"

java -version
echo "JAVA_HOME=${JAVA_HOME}"

#----------------------------------------
# 5. Rutas del proyecto
#----------------------------------------
project="quality-assessment"
home="/home/idalia/transcriptomics-ReMBE/dia_1"
data="/scratch/taller/datos_oyamel"
scripts="/home/idalia/transcriptomics-ReMBE/dia_1/scripts"

ADAPTERS="${SHARED_ENV}/share/trimmomatic-0.41-0/adapters/TruSeq3-PE-2.fa"

if [ ! -f "$ADAPTERS" ]; then
    echo "[ERROR] No se encontró el archivo de adaptadores: $ADAPTERS" >&2
    echo "[INFO] Búscalo con: find ${SHARED_ENV} -name 'TruSeq3-PE-2.fa'" >&2
    exit 1
fi

#----------------------------------------
# 6. Estructura de salida
#----------------------------------------
mkdir -p "${home}/${project}"/{work,results,logs}
cd "${home}/${project}"

#----------------------------------------
# 7. Ejecutar Nextflow
#----------------------------------------
export NXF_WORK="${home}/${project}/work"
export NXF_LOG_FILE="${home}/${project}/logs/nextflow.log"

srun nextflow run "${scripts}/01_main.nf" \
    -profile cluster \
    --reads "${data}/RAW/*.fastq.gz" \
    --outdir "${home}/${project}/results" \
    --adapters "${ADAPTERS}" \
    -with-report "${home}/${project}/logs/report.html" \
    -with-trace  "${home}/${project}/logs/trace.txt" \
    -with-timeline "${home}/${project}/logs/timeline.html" \
    -resume

echo "[INFO] Quality assessment completado."
echo "[INFO] Resultados en: ${home}/${project}/results"
echo "[INFO] Informe MultiQC: ${home}/${project}/results/multiqc/multiqc_report.html"