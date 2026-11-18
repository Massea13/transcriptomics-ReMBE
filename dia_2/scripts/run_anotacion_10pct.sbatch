#!/bin/bash
#SBATCH --job-name=anot_10pct
#SBATCH --output=anot_10pct_%j.out
#SBATCH --error=anot_10pct_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --partition=short

set -e

cd "$SLURM_SUBMIT_DIR"

#Anotación con TransDecoder y eggnog_env

# 1. Crear ligas simbólicas al ensamblado previo
mkdir -p ensamblado_10pct anotacion_10pct
ln -sf /scratch/marco/fastq/bioproject_985763/Trinity_10pct_final.fasta ensamblado_10pct/
ln -sf /scratch/marco/fastq/bioproject_985763/Trinity_10pct_final.gene_trans_map ensamblado_10pct/

cd anotacion_10pct

FASTA="../ensamblado_10pct/Trinity_10pct_final.fasta"
GENE_MAP="../ensamblado_10pct/Trinity_10pct_final.gene_trans_map"

# 2. Configurar entorno trinity_taller y localizar el binario real de TransDecoder
TD_ENV="/opt/ohpc/pub/apps/spack/local/linux-x86_64/miniconda3-24.7.1-g6mq4vwqwagz26lgp7shmzv6zyfrc5ui/envs/trinity_taller"
export PATH="$TD_ENV/bin:$PATH"

TD_LONG=$(find -L "$TD_ENV" -type f -name "TransDecoder.LongOrfs" 2>/dev/null | head -n 1)
TD_PRED=$(find -L "$TD_ENV" -type f -name "TransDecoder.Predict" 2>/dev/null | head -n 1)

if [ -z "$TD_LONG" ] || [ -z "$TD_PRED" ]; then
    echo "ERROR: No se encontró un ejecutable válido de TransDecoder." >&2
    exit 1
fi

echo "Usando TransDecoder: $TD_LONG"

echo "=== Fase 1: TransDecoder.LongOrfs ==="
"$TD_LONG" -t "$FASTA" --gene_trans_map "$GENE_MAP"

echo "=== Fase 2: TransDecoder.Predict ==="
"$TD_PRED" -t "$FASTA" --single_best_only

# 3. Activar eggnog_env para anotación funcional
export PATH="/scratch/marco/envs/eggnog_env/bin:$PATH"

TMP_EGG="/tmp/${USER}_eggnog_${SLURM_JOB_ID}"
mkdir -p "$TMP_EGG"
trap 'rm -rf "$TMP_EGG"' EXIT

PEP_FILE="Trinity_10pct_final.fasta.transdecoder.pep"

echo "=== Fase 3: eggNOG-mapper ==="
emapper.py -i "$PEP_FILE" \
           --itype proteins \
           -m diamond \
           --data_dir /scratch/marco/DIAMOND \
           --cpu $SLURM_CPUS_PER_TASK \
           --dbmem \
           --scratch_dir "$TMP_EGG" \
           --temp_dir "$TMP_EGG" \
           -o Trinity_10pct_funcional \
           --override

echo "Anotación estructural y funcional finalizada."
