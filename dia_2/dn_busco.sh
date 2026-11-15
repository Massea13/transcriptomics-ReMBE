#!/bin/bash --login
#SBATCH -J busco
#SBATCH -c 4   # (or --ntasks=) where numcores is between 2 and 168.
#SBATCH -o busco_plot%x.%j.out        # Standard output log
#SBATCH -e busco_plot%x.%j.err        # Standard error log
#SBATCH -t 1-0        # Wallclock limit (days-hours). Required!
                      # Max permitted is 7 days (7-0).


# Clean env
module purge
module load spack
module load miniconda3
source activate trinity_taller

FASTA_FILE="/scratch/taller/datos_pez/transcriptoma_completo/GKNI01.1.fsa_nt.gz"
OUT_DIR="calidad_ensamble"

# Comandos de BUSCO (se hace una sola vez, tarda ~ 1 hora)
busco -f -c 4 \
       -o ${OUT_DIR} \
       -i ${FASTA_FILE} \
       -l vertebrata_odb10 \
       -m transcriptome

# Comandos para generar los plots de BUSCO (se hace una sola vez, tarda ~ 1 minuto)
mkdir -p busco_plots
cp ${OUT_DIR}/short_summary*.json busco_plots/
busco --plot busco_plots --plot_percentages

