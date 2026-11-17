#!/bin/bash
#SBATCH --job-name=qc_10pct
#SBATCH --output=qc_10pct_%j.out
#SBATCH --error=qc_10pct_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=short

set -e

cd "$SLURM_SUBMIT_DIR"
mkdir -p qc_resultados

FILES=$(find -L . -maxdepth 2 -type f -name "*.fastq.gz" | sort)

if [ -z "$FILES" ]; then
    echo "ERROR: No se encontraron archivos .fastq.gz." >&2
    exit 1
fi

unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL
source /opt/ohpc/pub/apps/spack/local/linux-x86_64/miniconda3-24.7.1-g6mq4vwqwagz26lgp7shmzv6zyfrc5ui/etc/profile.d/conda.sh
conda activate trinity_taller
fastqc -t $SLURM_CPUS_PER_TASK -o qc_resultados $FILES
conda deactivate

conda activate rembe_qc_taller
multiqc -f qc_resultados -o qc_resultados/multiqc_report

echo "Control de calidad y reporte MultiQC finalizados."
