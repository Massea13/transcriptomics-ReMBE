- ## Sesión 1: Fundamentos, Entorno de Trabajo y Control de Calidad (4 horas)

  **Objetivo**: Establecer las bases conceptuales y técnicas del análisis de RNA-seq, familiarizarse con el entorno HPC/Nextflow y realizar el control de calidad de los datos crudos.

  **Teoría (1.5 horas)**:

#### 1. Introducción al análisis de RNA-seq:

  - Diseño experimental.
  - Contrastes
  - Réplicas biológicas y técnicas.
  - Efecto del lote de secuenciación (Batch effect)
  - Secuenciación Single-end vs. Paired-end, direccionalidad (strandedness)
  - Formatos de archivo: FASTQ (estructura y calidad), FASTA (secuencias de referencia), GTF/GFF (anotación genómica).

    

#### 2. Estrategias de análisis comparadas

  - Cuantificación por referencia: Mapeo de lecturas a un genoma o transcriptoma de referencia conocido.

  - Cuantificación sin referencia (*De novo*): Ensamblaje de lecturas para construir un transcriptoma *de novo*, para organismos sin genoma de referencia o no modelo.

##### 2.1 Presentación de los conjuntos de datos**

  - *Abies religiosa* (Oyameles)
  - *Genypterus maculatus* (Pex pistola)

   

#### 3. Flujo de trabajo general de un pipeline de RNA-seq.**

  ***Introducción a Conda, Nextflow y HPC:***

##### 3.1 ¿Qué es Conda?
Conda es un **gestor de paquetes y de entornos (package and environment manager)** de código abierto. Se usa principalmente en bioinformática, ciencia de datos y computación científica para instalar software y sus dependencias de forma aislada y reproducible [source](https://www.anaconda.com/about-us#).

**Gestor de paquetes**
Instala software ya compilado ("binarios") desde repositorios como **Anaconda**, **conda-forge** o **bioconda**, resolviendo automáticamente las dependencias (bibliotecas, versiones compatibles, etc.). Ejemplo

bash

```
conda install -c bioconda fastqc multiqc
```

**Gestor de entornos**
Crea espacios aislados donde cada proyecto tiene sus propias versiones de software, sin interferir entre sí ni con el sistema.

bash
```
conda create -n rembe_qc -c bioconda -c conda-forge fastqc multiqc
conda activate rembe_qc    # entras al entorno
conda deactivate           # sales del entorno
```

##### 3.2 ¿Qué es Nextflow?
Nextflow es un motor de orquestación de flujos de trabajo y un lenguaje específico de dominio (DSL: Domain Specific Languge) de código abierto, diseñado principalmente para facilitar la escritura y ejecución de flujos de trabajo computacionales intensivos en datos . En el ámbito de la bioinformática, **se utiliza comúnmente para conectar diversas herramientas de análisis (como FastQC, MultiQC, alineadores, etc.) ** y construir procesos de análisis reproducibles, portables y escalables



**1. Si quisieras hornear un pastel, ¿qué tareas ejecutaría Conda y cuáles corresponderían a Nextflow?**

<!--
Lista de ingredientes exactos → environment.yml  (parte de Conda)
Almacén que provee los ingredientes       → 		Conda
Receta (pasos y orden) + chef       → 			Nextflow
-->

##### Conda vs Nextflow

- **Conda** → gestiona **software** (qué programas y qué versiones hay instalados).
- **Nextflow** → gestiona **ejecución** (en qué orden se ejecutan los programas y en qué máquina).


| Aspecto                 | **Conda**                                                    | **Nextflow**                                                 |
| :---------------------- | :----------------------------------------------------------- | :----------------------------------------------------------- |
| **Qué es**              | Gestor de paquetes y entornos                                | Motor de orquestación de flujos de trabajo                   |
| **Qué gestiona**        | **Software**: instala/aisla versiones de programas (FastQC, Python, Java…) | **Ejecución**: orden, paralelismo, dependencias entre pasos, reintentos |
| **Unidad de trabajo**   | Un **paquete** (p. ej. `fastqc=0.12.1`)                      | Un **proceso** (p. ej. ejecutar FastQC sobre una muestra)    |
| **Aísla**               | Entornos de software (`rembe_qc`, `rembe_align`…)            | Pasos del pipeline y su ejecución (local, SLURM, nube)       |
| **Lenguaje**            | Comandos de shell (`conda install/create/activate`)          | DSL propio (`main.nf`: `process`, `workflow`, `channel`)     |
| **Escala**              | Un entorno por proyecto/máquina                              | Miles de tareas en paralelo en un clúster o nube             |
| **Reproducibilidad**    | Fija **versiones de software** (`environment.yml`)           | Fija **el flujo completo** + puede fijar software vía contenedores/Conda |
| **Reanudar tras fallo** | ❌ No aplica                                                  | ✅ `-resume` (omite pasos ya hechos)                          |
| **Se ejecuta en**       | Tu shell / nodo                                              | Local, SLURM, SGE, Kubernetes, AWS…                          |



#### 4. ¿Qué es la computación de alto rendimiento (HPC)?

La HPC es una tecnología que utiliza clústeres de potentes procesadores que trabajan en paralelo para procesar conjuntos de datos masivos y multidimensionales y resolver problemas complejos a velocidades extremadamente altas.

La HPC resuelve en tiempo real algunos de los problemas informáticos más complejos de la actualidad. Los sistemas HPC suelen funcionar a velocidades más de un millón de veces superiores a las de los sistemas básicos más rápidos de ordenadores de sobremesa, portátiles o servidores [fuente](https://www.ibm.com/es-es/think/topics/hpc).



**1. Pregunta.  Cual es la ventaja de usar  conda y/o nextflow en un entono  como los HPC?**

<!--
Conda: Resuelve el problema de las dependencias de software

La ventaja de Conda en HPC es que **proporciona entornos de software aislados y reproducibles**, permitiéndote evitar el sistema de "módulos" (modules) del clúster, que es lento y a veces restrictivo.

Nextflow: Resuelve el problema de la orquestación de flujos de trabajo

La ventaja de Nextflow en HPC es que **abstrae la complejidad de la programación de trabajos**, permitiéndote gestionar la paralización, el manejo de errores y la migración entre plataformas con código simple.
--> 



#### 5. Conceptos de paralelización, reproducibilidad y estructura básica de un pipeline (main.nf, nextflow.config).







  **Práctica en HPC con Nextflow (2.5 horas):**

  ​    1. Configuración del entorno: Acceso al cluster, carga de módulos (Nextflow, Java) y examinación del conjunto de datos (FASTQ, genoma de referencia y anotación) 

  ​    2. Control de Calidad Inicial: Ejecución de FastQC sobre los archivos FASTQ crudos para evaluar la calidad de las bases, el contenido de GC, la presencia de adaptadores y las secuencias sobrerrepresentadas.

  ​    3. Limpieza y Poda de Lecturas: Uso de Trimmomatic 

  ​      Eliminar adaptadores.

  ​      Recortar bases de baja calidad (umbral Q<30).

  ​      Filtrar lecturas demasiado cortas.

  ​      Manejo correcto de datos paired-end

  ​    4. Control de Calidad Post-Limpieza: Reejecución de FastQC en los archivos filtrados para verificar la mejora de la calidad.

      5. Reporte Agregado: Uso de MultiQC para compilar todos los reportes de FastQC en un solo informe HTML interactivo.

  












**Referencias**
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12795997/pdf/main.pdf

- https://monashbioinformaticsplatform.github.io/RNAseq_workshop_2024/experimental-design.html
- https://researchgold.org/blog/rna-sequencing-analysis
- https://hbctraining.github.io/Intro-to-rnaseq-fasrc-salmon-flipped/lessons/02_experimental_planning_considerations.html
- [IBM](https://www.ibm.com/es-es/think/topics/hpc)