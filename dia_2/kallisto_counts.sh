#!/bin/bash --login
#SBATCH -p multicore
#SBATCH -J kallisto_counts
#SBATCH --array=1-28
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=kallisto_%A_%a.out
#SBATCH --error=kallisto_%A_%a.err


module load apps/bioinf
module load apps/kallisto/0.43.0

# -------------------------
# Paths
# -------------------------
INDEX="trinity_denovo/Trinity_pez_dn.idx"
INPUT_DIR="trinity_denovo"
OUTPUT_DIR="kallisto_results"

THREADS=$SLURM_CPUS_PER_TASK

# -------------------------
# Setup
# -------------------------
mkdir -p "$OUTPUT_DIR"

# -------------------------
# Select sample
# -------------------------
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

READ1="${INPUT_DIR}/${SAMPLE}_clean_R1.fastq.gz"
READ2="${INPUT_DIR}/${SAMPLE}_clean_R2.fastq.gz"
OUTDIR="${OUTPUT_DIR}/${SAMPLE}"

# -------------------------
# Indexing just once
# -------------------------

kallisto index -i Trinity_pez_dn.idx Trinity-DN.fasta

# -------------------------
# Run kallisto
# -------------------------
kallisto quant \
  -i "$INDEX" \
  -o "$OUTDIR" \
  -t "$THREADS" \
  "$READ1" "$READ2"