#!/bin/bash
#SBATCH -J QC                 # Nombre del job
#SBATCH -p short              # Partición
#SBATCH -N 1                  # Nodos
#SBATCH -c 4                  # Núcleos (subido de 1 a 4 para FastQC/MultiQC)
#SBATCH --mem=8G              # Memoria (añadido, buena práctica)
#SBATCH -t 02:00:00           # Tiempo máximo (añadido)
#SBATCH -o QC-%j.out
#SBATCH -e QC-%j.err

set -euo pipefail   # Detiene el script ante cualquier error

#----------------------------------------
# 1. Cargar módulos del cluster
#----------------------------------------
module load spack/1.0.2
module load jdk/21.0.2-gcc-12.2.0-linux-almalinux9-x86_64-tsyhzjs
module load nextflow

# Forzar el uso del Java del módulo del clúster, sobrescribiendo JAVA_HOME de Conda
export JAVA_HOME=$(dirname $(dirname $(readlink -f $(which java))))
export JAVA_CMD="${JAVA_HOME}/bin/java"
export PATH="${JAVA_HOME}/bin:${PATH}"

# Verificar
java -version
echo "JAVA_HOME=${JAVA_HOME}"

#----------------------------------------
# 2. Instalar Miniconda si no existe
#----------------------------------------
CONDA_DIR="$HOME/miniconda3"

if [ ! -d "$CONDA_DIR" ]; then
    echo "[INFO] Miniconda no encontrada. Instalando en $CONDA_DIR ..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$CONDA_DIR"
    rm /tmp/miniconda.sh
    # Inicializar conda en bashrc (una sola vez)
    "$CONDA_DIR/bin/conda" init bash
fi

# Cargar conda en esta sesión
eval "$("$CONDA_DIR/bin/conda" shell.bash hook)"

#----------------------------------------
# 3. Crear el ambiente rembe_qc si no existe e instalar herramientas
#----------------------------------------
if ! conda env list | grep -qE "^rembe_qc\s"; then
    echo "[INFO] Creando ambiente conda 'rembe_qc' ..."
    conda create -y -n rembe_qc -c bioconda -c conda-forge fastqc multiqc
else
    echo "[INFO] El ambiente 'rembe_qc' ya existe. Verificando paquetes ..."
    conda install -y -n rembe_qc -c bioconda -c conda-forge fastqc multiqc
fi

conda activate rembe_qc

# Verificación rápida
fastqc --version
multiqc --version

#----------------------------------------
# 4. Organización del proyecto
#----------------------------------------
project="quality-assessment"
home="/home/idalia/transcriptomics-ReMBE/dia_1"
data="/home/idalia/rembe_common/data"
scripts="{$home}/scripts"

# Crear estructura de trabajo
mkdir -p "${home}/${project}"/{work,results,logs}
cd "${home}/${project}"

#----------------------------------------
# 5. Ejecutar Nextflow (dentro del job, con srun para que use los cores)
#----------------------------------------
export NXF_WORK="${scripts}/${project}/work"
export NXF_LOG_FILE="${home}/${project}/logs/nextflow.log"

srun nextflow run ./main.nf \
    -profile cluster \
    --reads  "${data}/*.fastq.gz" \
    --outdir "${home}/${project}/results" \
    -with-report "${home}/${project}/logs/report.html" \
    -with-trace  "${home}/${project}/logs/trace.txt" \
    -with-timeline "${home}/${project}/logs/timeline.html" \
    -resume

echo "[INFO] Quality assessment completado."




#Remover  todos los módulos cargados

#conda deactivate
#module purge 