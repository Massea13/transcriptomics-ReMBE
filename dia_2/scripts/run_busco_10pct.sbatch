#!/bin/bash
#SBATCH --job-name=busco_10pct
#SBATCH --output=busco_10pct_%j.out
#SBATCH --error=busco_10pct_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --partition=short

set -e

cd "$SLURM_SUBMIT_DIR"

unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL
source /opt/ohpc/pub/apps/spack/local/linux-x86_64/miniconda3-24.7.1-g6mq4vwqwagz26lgp7shmzv6zyfrc5ui/etc/profile.d/conda.sh
conda activate trinity_taller
export PATH="/opt/ohpc/pub/apps/spack/local/linux-x86_64/miniconda3-24.7.1-g6mq4vwqwagz26lgp7shmzv6zyfrc5ui/envs/trinity_taller/bin:$PATH"

# Buscar el FASTA en ensamblado_10pct/ o en el directorio actual
if [ -s "ensamblado_10pct/Trinity_10pct_final.fasta" ]; then
    FASTA=$(readlink -f ensamblado_10pct/Trinity_10pct_final.fasta)
elif [ -s "Trinity_10pct_final.fasta" ]; then
    FASTA=$(readlink -f Trinity_10pct_final.fasta)
else
    echo "ERROR: No se encontró Trinity_10pct_final.fasta en ensamblado_10pct/ ni en $SLURM_SUBMIT_DIR" >&2
    exit 1
fi

BUSCO_DB="/scratch/marco/busco_downloads"
export TMPDIR="/tmp/${USER}_busco_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

echo "Evaluando ensamblado: $FASTA"

busco -i "$FASTA" \
      -o busco_10pct \
      --out_path "$SLURM_SUBMIT_DIR" \
      -m transcriptome \
      -l metazoa_odb10 \
      --offline \
      --cpu $SLURM_CPUS_PER_TASK \
      --download_path "$BUSCO_DB" \
      -f

echo "Evaluación con BUSCO finalizada."
