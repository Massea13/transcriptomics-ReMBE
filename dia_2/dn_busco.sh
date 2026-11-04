#!/bin/bash --login
#SBATCH -p multicore
#SBATCH -J busco
#SBATCH -n 4   
#SBATCH -o busco_dn_%x.%j.out        # Standard output log
#SBATCH -e busco_dn_%x.%j.err        # Standard error log
#SBATCH -t 1-0        

#/mnt/iusers01/fatpou01/bmh01/m23015df/trimmed/cleaned_kraken_2/filtered

# Load any required modulefiles. A purge is used to start with a clean environment.
#module load apps/gcc/angsd/936

# Clean env
#conda deactivate
module purge
module load miniforge/24.11.2

conda activate busco_environment

# Paths to the input files and output directories
INPUT_FILE="trinity_denovo/Trinity-DN.fasta"
OUTPUT_DIR="trinity_denovo/busco_pez_results"

busco -f -c 16 -o $OUTPUT_DIR -i $INPUT_FILE -l vertebrata_odb10 -m transcriptome

generate_plot.py -wd $OUTPUT_DIR
