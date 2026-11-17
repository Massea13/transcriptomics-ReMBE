#!/bin/bash
#SBATCH --job-name=trim_10pct
#SBATCH --output=trim_10pct_%j.out
#SBATCH --error=trim_10pct_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --partition=short

set -e

source /opt/ohpc/pub/apps/spack/local/linux-x86_64/miniconda3-24.7.1-g6mq4vwqwagz26lgp7shmzv6zyfrc5ui/etc/profile.d/conda.sh
conda activate trinity_taller

cd "$SLURM_SUBMIT_DIR"

# Reservar 24 GB de heap explícitos para Java (evita OutOfMemoryError en GzipParallelCompressor)
export _JAVA_OPTIONS="-Xmx24G"

INPUT_DIR="subsampled_10pct"
mkdir -p trimmed_reads/paired trimmed_reads/unpaired trimmed_reads/logs

# Localizar adaptadores Illumina dentro del entorno
ADAPTERS=$(find "$CONDA_PREFIX" -name "TruSeq3-PE-2.fa" -o -name "TruSeq3-PE.fa" | head -n 1)

if [ -z "$ADAPTERS" ]; then
    echo "ERROR: No se encontró el archivo de adaptadores TruSeq3-PE en $CONDA_PREFIX" >&2
    exit 1
fi

echo "Usando adaptadores: $ADAPTERS"

for r1 in ${INPUT_DIR}/*_1.fastq.gz; do
    base=$(basename "$r1" _1.fastq.gz)
    r2="${INPUT_DIR}/${base}_2.fastq.gz"

    echo "Procesando muestra: $base..."

    trimmomatic PE -threads $SLURM_CPUS_PER_TASK -phred33 \
        "$r1" "$r2" \
        "trimmed_reads/paired/${base}_1.paired.fastq.gz" \
        "trimmed_reads/unpaired/${base}_1.unpaired.fastq.gz" \
        "trimmed_reads/paired/${base}_2.paired.fastq.gz" \
        "trimmed_reads/unpaired/${base}_2.unpaired.fastq.gz" \
        ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36 \
        2> "trimmed_reads/logs/${base}_trimmomatic.log"
done

echo "--- Resumen de retención de lecturas ---"
grep "Input Read Pairs" trimmed_reads/logs/*_trimmomatic.log
echo "Limpieza con Trimmomatic finalizada exitosamente."
