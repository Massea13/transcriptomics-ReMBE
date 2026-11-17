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
# 1. Cargar módulos del cluster
#----------------------------------------
module load spack
module load miniconda3
module load jdk/21.0.2-gcc-12.2.0-linux-almalinux9-x86_64-tsyhzjs
module load nextflow

# Guardar el Java del módulo (21) ANTES de activar conda,
# porque el entorno trae su propio Java 25, que Nextflow no acepta
MODULE_JAVA_HOME=$(dirname "$(dirname "$(readlink -f "$(which java)")")")
echo "[INFO] Java del módulo: ${MODULE_JAVA_HOME}"

#----------------------------------------
# 2. Activar el entorno COMPARTIDO por ruta absoluta
#----------------------------------------
SHARED_ENV="/opt/ohpc/pub/apps/spack/local/linux-x86_64/miniconda3-24.7.1-g6mq4vwqwagz26lgp7shmzv6zyfrc5ui/envs/rembe_qc"

if [ ! -d "$SHARED_ENV" ]; then
    echo "[ERROR] El entorno compartido no existe: $SHARED_ENV" >&2
    exit 1
fi

# Los scripts de activación de conda pueden usar variables sin definir
set +u
source activate "$SHARED_ENV"
set -u

# Verificar que las herramientas están disponibles
fastqc --version
multiqc --version
trimmomatic -version

#----------------------------------------
# 3. Forzar Java del módulo (sobrescribe JAVA_HOME del entorno compartido)
#----------------------------------------
export JAVA_HOME="${MODULE_JAVA_HOME}"
export JAVA_CMD="${JAVA_HOME}/bin/java"
export NXF_JAVA_HOME="${JAVA_HOME}"
export PATH="${JAVA_HOME}/bin:${PATH}"

java -version
echo "JAVA_HOME=${JAVA_HOME}"

#----------------------------------------
# 4. Rutas del proyecto
#----------------------------------------
project="quality-assessment"
# Se adapta solo a cada usuario (daniela, taller18, ...)
home="${HOME}/transcriptomics-ReMBE/dia_1"
data="/scratch/taller/datos_oyamel"
scripts="${home}/scripts"

# Busca los adaptadores sin depender de la versión de Trimmomatic
ADAPTERS=$(find "${SHARED_ENV}/share" -name 'TruSeq3-PE-2.fa' 2>/dev/null | head -1 || true)

if [ -z "$ADAPTERS" ] || [ ! -f "$ADAPTERS" ]; then
    echo "[ERROR] No se encontró TruSeq3-PE-2.fa dentro de ${SHARED_ENV}/share" >&2
    exit 1
fi
echo "[INFO] Adaptadores: $ADAPTERS"

#----------------------------------------
# 5. Estructura de salida
#----------------------------------------
mkdir -p "${home}/${project}"/{work,results,logs}
cd "${home}/${project}"

#----------------------------------------
# 6. Ejecutar Nextflow
#----------------------------------------
export NXF_WORK="${home}/${project}/work"
export NXF_LOG_FILE="${home}/${project}/logs/nextflow.log"

# Config extra que obliga a TODOS los procesos a usar el entorno compartido,
# sin importar qué ruta de conda traiga escrita 01_main.nf
# (withName tiene prioridad sobre la directiva `conda` del script)
OVERRIDE_CONFIG="${home}/${project}/logs/conda_override.config"
cat > "$OVERRIDE_CONFIG" <<EOF
conda.enabled = true
process {
    withName: '.*' {
        conda = '${SHARED_ENV}'
    }
}
EOF

srun nextflow run "${scripts}/01_main.nf" \
    -profile cluster \
    -c "$OVERRIDE_CONFIG" \
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
