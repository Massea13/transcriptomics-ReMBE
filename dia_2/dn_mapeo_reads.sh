#!/bin/bash --login
#SBATCH -p multicore
#SBATCH -J star_array
#SBATCH --cpus-per-task=8
#SBATCH --mem=90G
#SBATCH -t 1-0
#SBATCH --array=1-28%4 # have to fix the number of samples in samples.txt
#SBATCH -o star_%A_%a.out
#SBATCH -e star_%A_%a.err


module load apps/bioinf
module load apps/binapps/star/2.7.2b

# -------------------------
# Paths
# -------------------------
INPUT_DIR="reads_trimmed/"
OUTPUT_DIR="star_results_reads_mapeados"
THREADS=8

# -------------------------
# Setup
# -------------------------
mkdir -p "$OUTPUT_DIR"

# -------------------------
# Select sample for this array task
# -------------------------
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

READ1="${INPUT_DIR}/${SAMPLE}_clean_R1.fastq.gz"
READ2="${INPUT_DIR}/${SAMPLE}_clean_R2.fastq.gz"
PREFIX="${OUTPUT_DIR}/${SAMPLE}_"

# -------------------------
# Run STAR
# -------------------------
echo "[$(date)] START STAR for $SAMPLE"
echo "Threads: $THREADS"

STAR \
  --genomeDir "$GENOME_DIR" \
  --readFilesIn "$READ1" "$READ2" \
  --readFilesCommand zcat \
  --runThreadN "$THREADS" \
  --outFileNamePrefix "$PREFIX" \
  --outSAMtype BAM SortedByCoordinate

echo "[$(date)] DONE STAR for $SAMPLE"