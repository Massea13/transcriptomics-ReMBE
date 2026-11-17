#!/bin/bash
#SBATCH --job-name=trinity_trim
#SBATCH --output=trinity_trim_%j.out
#SBATCH --error=trinity_trim_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --partition=long

set -e

source /opt/ohpc/pub/apps/spack/local/linux-x86_64/miniconda3-24.7.1-g6mq4vwqwagz26lgp7shmzv6zyfrc5ui/etc/profile.d/conda.sh
conda activate trinity_taller

cd "$SLURM_SUBMIT_DIR"

# Construir listas ordenadas separadas por coma de R1 y R2 limpios
LEFT_READS=$(ls trimmed_reads/paired/*_1.paired.fastq.gz | sort | paste -sd "," -)
RIGHT_READS=$(ls trimmed_reads/paired/*_2.paired.fastq.gz | sort | paste -sd "," -)

if [ -z "$LEFT_READS" ] || [ -z "$RIGHT_READS" ]; then
    echo "ERROR: No se encontraron lecturas pareadas en trimmed_reads/paired/" >&2
    exit 1
fi

# Directorio temporal local en el nodo para evitar cuellos de botella de red en /scratch
TMP_WORK="/tmp/${USER}_trinity_trim_${SLURM_JOB_ID}"
mkdir -p "$TMP_WORK"

# Limpiar /tmp automáticamente al terminar o cancelar el trabajo
trap 'rm -rf "$TMP_WORK"' EXIT

echo "Iniciando ensamblaje con Trinity en $TMP_WORK..."
echo "Muestras Left: $LEFT_READS"
echo "Muestras Right: $RIGHT_READS"

Trinity --seqType fq \
        --max_memory 50G \
        --CPU $SLURM_CPUS_PER_TASK \
        --left "$LEFT_READS" \
        --right "$RIGHT_READS" \
        --min_kmer_cov 2 \
        --no_salmon \
        --full_cleanup \
        --output "$TMP_WORK/trinity_out"

# Rescatar el ensamblado final y el mapa gen-transcrito hacia /scratch
FASTA_OUT=$(find "$TMP_WORK" -maxdepth 2 -type f -name "*.Trinity.fasta" | head -n 1)
MAP_OUT=$(find "$TMP_WORK" -maxdepth 2 -type f -name "*.gene_trans_map" | head -n 1)

if [ -z "$FASTA_OUT" ] || [ ! -s "$FASTA_OUT" ]; then
    echo "ERROR: No se encontró el archivo .Trinity.fasta generado en $TMP_WORK" >&2
    exit 1
fi

cp "$FASTA_OUT" "$SLURM_SUBMIT_DIR/Trinity_trimmed_final.fasta"
[ -n "$MAP_OUT" ] && cp "$MAP_OUT" "$SLURM_SUBMIT_DIR/Trinity_trimmed_final.gene_trans_map"

# Generar estadísticas básicas del ensamblado si TrinityStats.pl está disponible
if command -v TrinityStats.pl &> /dev/null; then
    TrinityStats.pl "$SLURM_SUBMIT_DIR/Trinity_trimmed_final.fasta" > "$SLURM_SUBMIT_DIR/Trinity_trimmed_stats.txt"
fi

echo "Ensamblaje finalizado y guardado en: $SLURM_SUBMIT_DIR/Trinity_trimmed_final.fasta"
