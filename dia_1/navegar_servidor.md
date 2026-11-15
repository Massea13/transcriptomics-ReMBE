## Guía rápida para navegar en el servidor
#### 1. ¿Qué es SLURM?
SLURM (Simple Linux Utility for Resource Management) es un sistema de gestión de trabajos y recursos para clústeres Linux. Permite encolar y ejecutar trabajos (jobs) en nodos de cómputo, asignando recursos (CPU, memoria, GPU, tiempo) de forma eficiente y compartida entre múltiples usuarios.

## Comandos básicos de navegación

### Navegación y exploración de archivos

| Comando           | Descripción                                                  |
| :---------------- | :----------------------------------------------------------- |
| `ls`              | Lista los archivos y directorios del directorio actual.      |
| `ls -lth`         | Lista en formato largo (`-l`), ordenado por fecha de modificación (`-t`), con tamaños legibles (`-h`). Ideal para ver qué se modificó recientemente. |
| `cd <directorio>` | Cambia al directorio especificado.                           |
| `cd ..`           | Sube **un directorio arriba** (directorio padre).            |
| `cd .`            | Referencia al **directorio actual** (poco usado para moverse, pero útil en rutas). |
| `cat <archivo>`   | Muestra todo el contenido de un archivo en pantalla.         |
| `less <archivo>`  | Abre el archivo en un visor paginado (navega con flechas, sale con `q`). Mejor para archivos grandes. |
| `head <archivo>`  | Muestra las primeras 10 líneas de un archivo (usá `head -n 20` para más líneas). |

### SLURM

| Comando              | Descripción                                                  |
| :------------------- | :----------------------------------------------------------- |
| `sbatch <script.sh>` | Envía un script de trabajo a la cola de SLURM para su ejecución. |

> 💡 Un script típico de `sbatch` incluye directivas como `#SBATCH --job-name`, `#SBATCH --cpus-per-task`, `#SBATCH --mem`, `#SBATCH --time`, seguidas de los comandos a ejecutar.

### Conda (entornos)

| Comando                           | Descripción                                                  |
| :-------------------------------- | :----------------------------------------------------------- |
| `conda`                           | Gestor de paquetes y entornos (verifica que esté disponible con `conda --version`). |
| `conda env list`                  | Lista todos los entornos Conda disponibles (también sirve `conda info --envs`). |
| `conda activate <nombre_entorno>` | Activa el entorno Conda especificado.                        |
| `conda deactivate`                | Desactiva el entorno actual y vuelve al base.                |