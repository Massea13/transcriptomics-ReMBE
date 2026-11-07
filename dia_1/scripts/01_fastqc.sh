#SBATCH -J QC			# Nombre
#SBATCH -p short		# Particion donde se ejecuta el job
#SBATCH -N 1			# Numero de nodos
#SBATCH -c 1                    # Numero de nucleos
#SBATCH -o QC-%j.out	# Nombre del archivo de salida
#SBATCH -e QC-%j.err	# Nombre del archivo de errores

# Cargar modulos necesarios 
module load spack/1.0.2
module load jdk/21.0.2-gcc-12.2.0-linux-almalinux9-x86_64-tsyhzjs
module load nextflow


# 1. Descargar el script de instalación de Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# 2. Ejecutar la instalación (se instalará en el directorio ~/miniconda3 por defecto)
bash Miniconda3-latest-Linux-x86_64.sh

#install fastqc and multiqc
conda create -n rembe_qc -c bioconda -c conda-forge fastqc multiqc









#Remover  todos los módulos cargados
#module purge 
