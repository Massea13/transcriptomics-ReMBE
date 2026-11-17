# Ensamblaje de novo con Trinity: calidad, anotación y cuantificación (teoría)

Este documento cubre la primera mitad de la teoría (45 min): cómo se ensambla un transcriptoma de novo con Trinity, cómo se evalúa y anota el ensamblaje, y cómo se cuantifica la expresión mapeando las lecturas contra él.

## Cómo encaja en las 1.5 h de teoría

| Parte | Tema | Tiempo | Dónde está |
| --- | --- | --- | --- |
| 1 | Ensamble de novo: introducción a Trinity | 15 min | Este documento, Tema 1 |
| 1 | Evaluación de la calidad y anotación del ensamblaje | 15 min | Este documento, Tema 2 |
| 1 | Cuantificación con un transcriptoma de novo: mapeo de lecturas contra el ensamblaje | 15 min | Este documento, Tema 3 |
| 2 | Consideraciones biológicas: cuándo usar de novo, limitaciones (isoformas, quimeras, redundancia) | 15 min | Documento "Cuantificación sin transcriptoma de referencia (ensamblaje de novo): teoría", Tema 1 |
| 2 | Pseudoalineamiento: Salmon y Kallisto | 15 min | Ese mismo documento, Tema 2 |
| 2 | Estructura de índices: qué archivos se necesitan y cómo se construyen | 15 min | Ese mismo documento, Tema 3 |

Los tres últimos temas de tu lista ya estaban en el documento anterior; aquí no los repito para no duplicar contenido. El orden de clase sugerido es este: primero el ensamblaje y su evaluación, después la cuantificación y, al final, las limitaciones y el detalle de Salmon, Kallisto e índices.

## Flujo general

```text
FASTQ crudo -> QC y recorte suave (FastQC, fastp)                          [Tema 1: entrada]
            -> ensamblaje de novo (Trinity)                                [Tema 1]
            -> Trinity.fasta + gene_trans_map
            -> evaluacion: N50, ExN50, BUSCO, back-mapping                 [Tema 2]
            -> anotacion: TransDecoder, DIAMOND/BLAST, Pfam, Trinotate     [Tema 2]
            -> mapeo/pseudoalineamiento de lecturas por muestra            [Tema 3]
            -> matrices de conteo y TPM (transcrito y gen)                 [Tema 3]
            -> expresion diferencial (DESeq2 / edgeR)
```

## Tema 1. Ensamble de novo: introducción a Trinity (15 min)

Trinity reconstruye transcritos completos a partir de lecturas Illumina sin genoma de referencia, usando grafos de De Bruijn y tres módulos que se aplican en secuencia: Inchworm, Chrysalis y Butterfly ([Grabherr et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3571712)).

**Reparto sugerido:** 5 min el problema y los tres módulos, 5 min entradas y parámetros, 5 min salidas, recursos y cómo correrlo.

### 1.1 El problema que resuelve

En RNA-seq no se cumple el supuesto de "un locus, un contig" ni el de cobertura uniforme: un ensamblador de transcriptomas debe producir un contig por transcrito (isoforma) y no por locus, con genes que se expresan miles de veces más que otros ([protocolo de Haas et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3875132)). Por eso se usa un ensamblador específico y no uno de genomas.

Idea de un **grafo de De Bruijn**: las lecturas se parten en k-mers; cada k-mer es un nodo y dos k-mers que se solapan en k−1 bases se unen; un transcrito es un camino en el grafo. En la práctica Trinity usa k = 25 ([artículo original](https://scispace.com/pdf/full-length-transcriptome-assembly-from-rna-seq-data-without-55o8f60wxt.pdf)).

### 1.2 Los tres módulos

| Módulo | Qué hace | Salida conceptual | Recurso que consume |
| --- | --- | --- | --- |
| Inchworm | Construye un diccionario de k-mers, descarta los k-mers que probablemente son errores y extiende de forma voraz los k-mers más frecuentes ([artículo original](https://scispace.com/pdf/full-length-transcriptome-assembly-from-rna-seq-data-without-55o8f60wxt.pdf)). Suele dar transcritos completos de la isoforma dominante y solo las porciones únicas de las isoformas alternativas ([CUNY HPC](https://wiki.csi.cuny.edu/cunyhpc/index.php/TRINITY)) | Contigs lineales | Memoria RAM |
| Chrysalis | Agrupa los contigs de Inchworm en clústeres y construye un grafo de De Bruijn completo por clúster; cada clúster representa la complejidad transcripcional de un gen o de un conjunto de genes que comparten secuencia ([wiki de Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Home)). Reparte las lecturas entre los grafos | Clústeres con su grafo | Memoria RAM |
| Butterfly | Procesa los grafos en paralelo, sigue los caminos que recorren lecturas y pares de lecturas, y reporta transcritos completos de isoformas alternativas, separando los que corresponden a palólogos ([MetaCentrum](https://wiki.metacentrum.cz/wiki/Trinity)) | Transcritos finales | Procesadores |

En la documentación de uso se agrupan en dos fases: la fase 1 (Jellyfish, Inchworm y Chrysalis) agrupa las lecturas, y la fase 2 (Butterfly) ensambla cada clúster ([guía de St Andrews](https://stab.st-andrews.ac.uk/wiki/index.php/Trinity)). Inchworm y Chrysalis dependen de la memoria y Butterfly del número de núcleos ([FASRC Harvard](https://docs.rc.fas.harvard.edu/kb/trinity/)), lo que importa al pedir recursos en el clúster.

### 1.3 Nombres de los contigs

```text
>TRINITY_DN1000_c0_g1_i2 len=1847 path=[...]
         |       |  |  |
         |       |  |  +-- isoforma dentro del gen
         |       |  +----- "gen" de Trinity (grupo de isoformas)
         |       +-------- sub-componente
         +---------------- componente (clúster de Chrysalis)
```

El "gen" de Trinity no es un gen anotado: es un grupo de transcritos con k-mers compartidos, y puede unir palólogos o separar alelos.

### 1.4 Entradas y parámetros importantes

| Parámetro | Para qué sirve | Sugerencia |
| --- | --- | --- |
| `--seqType fq` | Formato de las lecturas | FASTQ |
| `--samples_file` | Tabla de muestras: condición, réplica, R1 y R2, separadas por tabulador | Permite ensamblar varias muestras y condiciones juntas ([tutorial de UFS](https://docs.hpc.ufs.ac.za/training/transcriptomics_tutorial/transcriptomics_tutorial_1/)) |
| `--SS_lib_type` | Orientación de la librería estrandada (`RF` o `FR`) | Indicarla si la librería es estrandada; omitirla si no |
| `--CPU` y `--max_memory` | Recursos | El `--max_memory` debe ser menor que lo pedido en SLURM |
| `--min_contig_length` | Largo mínimo del contig | 200 por defecto |
| `--no_normalize_reads` | Desactiva la normalización in silico | Se deja activa con muchísimas lecturas; reduce memoria y tiempo |
| `--genome_guided_bam` | Ensamblaje guiado por un BAM alineado a un genoma | Cuando hay genoma, pero la anotación es pobre |
| `--output` | Directorio de salida | Su nombre debe contener la palabra `trinity` |
| `--full_cleanup` | Borra archivos intermedios | Ahorra disco, pero impide reanudar si el trabajo muere |

**Qué reads usar.** Se recomienda un recorte suave de calidad y de adaptadores (por ejemplo con fastp); un recorte agresivo descarta información y tiende a producir ensamblajes más fragmentados. Una guía basada en evidencia propone entre 20 y 40 millones de lecturas de un solo individuo cuando sea posible ([MacManes, preprint](https://www.biorxiv.org/content/10.1101/035642v1.full-text)).

**Una muestra o varias.** Para cuantificar y comparar condiciones, se ensambla una referencia única con lecturas de todas las condiciones (o de un subconjunto representativo) y se cuantifica cada muestra contra ella (Tema 3).

### 1.5 Salidas

| Archivo | Contenido |
| --- | --- |
| `<salida>.Trinity.fasta` | Todos los transcritos del ensamblaje |
| `<salida>.Trinity.fasta.gene_trans_map` | Tabla gen → transcrito, que se usa para cuantificación a nivel gen |
| `<salida>/` | Intermedios: k-mers, clústeres de Chrysalis, lecturas normalizadas |

### 1.6 Cómo correrlo

```bash
# samples.txt (tabuladores): condicion  replica  R1  R2
# control  ctrl_rep1  /ruta/ctrl1_R1.trim.fq.gz  /ruta/ctrl1_R2.trim.fq.gz
# trat     trat_rep1  /ruta/trat1_R1.trim.fq.gz  /ruta/trat1_R2.trim.fq.gz

# MUY PESADO: requiere sbatch en el clúster (p. ej. --mem=220G, -c 16, 72 h)
Trinity --seqType fq --samples_file samples.txt \
        --CPU 16 --max_memory 200G \
        --min_contig_length 200 \
        --output trinity_proyecto

# Revisar tamaño y conteo de contigs
grep -c '>' trinity_proyecto.Trinity.fasta
```

La orden `Trinity` no es interactiva ni debe correrse en el nodo de acceso. Si el trabajo se interrumpe, al relanzarlo con el mismo directorio de salida retoma desde el último punto, siempre que no se haya usado `--full_cleanup`.

### 1.7 Estado del proyecto y alternativas

Según su wiki, desde 2024 el proyecto Trinity ya no cuenta con financiamiento dedicado para desarrollo, mantenimiento o soporte ([wiki de Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Home)). Sigue siendo el estándar para este curso, pero conviene saberlo. Otros ensambladores de transcriptomas (de memoria, sin verificar aquí): rnaSPAdes, SOAPdenovo-Trans y Trans-ABySS; algunos flujos combinan varios ensambladores y filtran la redundancia con EvidentialGene ([Duarte et al., protocolo en Bio-protocol](https://bio-protocol.org/en/bpdetail?id=3912)).

## Tema 2. Evaluación de la calidad y anotación del ensamblaje (15 min)

Un ensamblaje de novo no tiene "respuesta correcta" con la que compararlo, así que su calidad se juzga con varias medidas complementarias (contigs, completitud de genes, representación de las lecturas) y después se le da significado biológico prediciendo proteínas y buscando homología.

**Reparto sugerido:** 8 min evaluación, 7 min anotación.

### 2.1 Medidas de calidad

El flujo habitual evalúa el mapeo de lecturas, los estadísticos Nx y ExN50, la completitud con BUSCO y, opcionalmente, métricas de errores de ensamblaje con QUAST ([tutorial del NYU Genomics Core](https://learn.gencore.bio.nyu.edu/denovo/transcriptome.html)).

| Medida | Qué responde | Cómo se obtiene | Qué esperar |
| --- | --- | --- | --- |
| Estadísticos básicos | ¿Cuántos transcritos y "genes" hay y qué largos son? | `TrinityStats.pl` ([tutorial OICR](https://bioinformaticsdotca.github.io/rnaseq_2017_tutorial6)) | Se reportan con todas las isoformas y con la isoforma más larga por gen |
| N50 | Longitud tal que la mitad de las bases está en contigs de ese largo o mayores | `TrinityStats.pl` | Poco informativo por sí solo en transcriptomas |
| ExN50 | N50 calculado solo con los transcritos que suman x % de la expresión total | `contig_ExN50_statistic.pl` sobre una matriz de expresión | Suele crecer hasta un máximo cerca de Ex90 y luego caer |
| BUSCO | ¿Están los genes de copia única esperados para el linaje, completos, fragmentados, ausentes o duplicados? | `busco -m transcriptome` | Completos altos; duplicados altos indican redundancia |
| Back-mapping | ¿Qué fracción de las lecturas se alinea al ensamblaje? | `bowtie2` o el log de Salmon | Valores altos; bajos indican ensamblaje incompleto o contaminación |
| Errores de ensamblaje | Quimeras, cobertura anormal, contigs mal ensamblados | TransRate o QUAST ([ejemplo de uso de TransRate y BUSCO](https://eel-pond.readthedocs.io/en/latest/4-evaluating-assembly.html)) | Revisar el estado de mantenimiento de la herramienta antes de adoptarla |

**BUSCO:** un gen se considera completo cuando su longitud queda dentro de dos desviaciones estándar de la longitud media de su grupo ([ejemplo de uso](https://eel-pond.readthedocs.io/en/latest/4-evaluating-assembly.html)).

**Por qué ExN50 y no solo N50:** un ensamblaje con muchos contigs cortos y poco expresados (ruido, fragmentos) puede reducir el N50 sin afectar a los genes relevantes, y uno redundante o quimérico puede inflarlo. El ExN50 pesa por expresión, por lo que refleja mejor lo que se va a cuantificar. Este detalle y los rangos típicos de la tabla están escritos de memoria.

```bash
# 1) Estadísticos del ensamblaje
TrinityStats.pl trinity_proyecto.Trinity.fasta > TrinityStats.txt

# 2) Completitud (elegir el linaje del grupo taxonómico)
busco -i trinity_proyecto.Trinity.fasta -l insecta_odb10 -m transcriptome -c 16 -o busco_asm

# 3) Back-mapping de lecturas al ensamblaje (por muestra)
bowtie2-build trinity_proyecto.Trinity.fasta asm_index
bowtie2 -p 8 -q --no-unal -k 20 -x asm_index -1 R1.trim.fq.gz -2 R2.trim.fq.gz 2> muestra.bowtie2.log
# la tasa global aparece como "overall alignment rate" en el log

# 4) ExN50 (después de cuantificar, Tema 3)
contig_ExN50_statistic.pl matriz.TMM.EXPR.matrix trinity_proyecto.Trinity.fasta > ExN50.stats
```

BUSCO y bowtie2 consumen tiempo y memoria: se lanzan con `sbatch`. `TrinityStats.pl` es ligero.

### 2.2 Anotación funcional

La anotación responde dos preguntas: qué transcritos codifican proteínas y a qué se parecen. El flujo más usado predice los marcos de lectura con TransDecoder y anota con Trinotate, que integra BLAST o DIAMOND contra Swiss-Prot, dominios Pfam con HMMER, EggNOG, TmHMM y SignalP ([tutorial del NYU Genomics Core](https://learn.gencore.bio.nyu.edu/denovo/transcriptome.html)).

| Paso | Herramienta | Qué produce | Costo |
| --- | --- | --- | --- |
| Marcos de lectura | TransDecoder (`LongOrfs` y `Predict`) | Proteínas predichas (`.pep`) y CDS | Moderado |
| Homología de proteínas | DIAMOND o BLAST (`blastp`) contra Swiss-Prot o una base de un grupo cercano | Mejor hit por proteína | Alto; `sbatch` |
| Homología de transcritos | DIAMOND o BLAST (`blastx`) | Mejor hit por transcrito, incluye los no codificantes | Alto; `sbatch` |
| Dominios | HMMER (`hmmscan`) contra Pfam | Dominios proteicos | Alto; `sbatch` |
| Integración | Trinotate | Tabla única con GO, KEGG, Pfam, EggNOG y más | Moderado |
| Enriquecimiento funcional | GOseq u otras | Términos GO enriquecidos entre los DEGs ([flujo de ejemplo](https://hackmd.io/@chenyanniii/Hy0a37H3_)) | Bajo |

TransDecoder se corre en dos pasos y el segundo usa la misma carpeta de salida que el primero; en el segundo se pueden aportar los resultados de homología y de dominios para retener marcos de lectura con esa evidencia ([ejemplo de uso](https://hackmd.io/@chenyanniii/Hy0a37H3_)).

```bash
# sbatch para ensamblajes completos
TransDecoder.LongOrfs -t trinity_proyecto.Trinity.fasta --gene_trans_map trinity_proyecto.Trinity.fasta.gene_trans_map

diamond blastp -d uniprot_sprot.dmnd -q trinity_proyecto.Trinity.fasta.transdecoder_dir/longest_orfs.pep \
  -o orfs_vs_sprot.outfmt6 --outfmt 6 --max-target-seqs 1 --evalue 1e-5 -p 16

hmmscan --cpu 16 --domtblout pfam.domtblout Pfam-A.hmm \
  trinity_proyecto.Trinity.fasta.transdecoder_dir/longest_orfs.pep > /dev/null

TransDecoder.Predict -t trinity_proyecto.Trinity.fasta \
  --retain_blastp_hits orfs_vs_sprot.outfmt6 --retain_pfam_hits pfam.domtblout

diamond blastx -d uniprot_sprot.dmnd -q trinity_proyecto.Trinity.fasta \
  -o transcripts_vs_sprot.outfmt6 --outfmt 6 --max-target-seqs 1 --evalue 1e-5 --more-sensitive -p 16
```

Trinotate tiene una instalación y una sintaxis que cambian entre versiones: se sigue la documentación de la versión instalada (Bryant et al. 2017 describe la herramienta; [referencia citada en el protocolo de Duarte](https://pmc.ncbi.nlm.nih.gov/articles/PMC7952955)).

**Qué transmitir sobre la anotación:**

- La mayoría de los transcritos de un organismo no modelo queda sin anotar; eso es normal y depende de la distancia a las especies con proteínas bien anotadas.
- La base de datos importa: Swiss-Prot es precisa pero pequeña, y una base de un grupo cercano da más hits.
- Anotar no valida el ensamblaje: un contig quimérico puede tener hits a dos proteínas distintas.
- Las etiquetas por homología son hipótesis, no hechos experimentales.

## Tema 3. Cuantificación con un transcriptoma de novo: mapeo de lecturas contra el ensamblaje (15 min)

Sin genoma ni anotación, el propio ensamblaje de Trinity es la referencia: las lecturas de cada muestra se mapean o pseudoalinean contra él, se estima cuántas vienen de cada transcrito y los resultados se agregan a "gen" de Trinity para compararlos entre condiciones.

**Reparto sugerido:** 5 min la lógica y las opciones de herramienta, 5 min el flujo paso a paso con comandos, 5 min matrices, normalización y errores comunes.

### 3.1 La lógica

1. Se construye **una sola referencia** (el ensamblaje) con lecturas de todas las condiciones.
2. Cada **muestra** se cuantifica por separado contra esa misma referencia, con sus lecturas recortadas.
3. Una lectura puede ser compatible con varios transcritos (isoformas, alelos, contigs redundantes); el estimador reparte esas lecturas de forma probabilística (algoritmo EM).
4. Por muestra se obtienen **conteos estimados**, **TPM** y **longitud efectiva** de cada transcrito.
5. Las salidas de todas las muestras se juntan en **matrices** de transcritos por muestras y, con el mapa gen-transcrito, de "genes" por muestras.
6. La matriz de conteos a nivel gen entra a DESeq2 o edgeR.

**Por qué la misma referencia para todas:** si cada muestra tuviera su propio ensamblaje, los "genes" no serían comparables entre muestras.

### 3.2 Opciones de herramienta

| Opción | Enfoque | Cuándo conviene |
| --- | --- | --- |
| RSEM con bowtie2 | Alinea las lecturas (bowtie2 permitiendo muchos alineamientos por lectura) y estima la abundancia con EM | Flujo clásico de Trinity; es el que usa el tutorial del NYU ([NYU Genomics Core](https://learn.gencore.bio.nyu.edu/denovo/transcriptome.html)); más lento y con archivos BAM grandes |
| Salmon | Mapeo rápido de k-mers con alineamiento selectivo y EM | Opción por defecto en este curso por velocidad y poco uso de disco |
| Kallisto | Pseudoalineamiento por grafo de De Bruijn del transcriptoma | Alternativa equivalente a Salmon; no detecta el tipo de librería |

Los detalles de Salmon y Kallisto (algoritmo, salidas, índices) están en el documento "Cuantificación sin transcriptoma de referencia (ensamblaje de novo): teoría". Aquí se presenta el flujo de extremo a extremo.

### 3.3 Flujo con los scripts de Trinity

Trinity incluye `align_and_estimate_abundance.pl`, que prepara la referencia y llama al método elegido (RSEM, kallisto o salmon), y `abundance_estimates_to_matrix.pl`, que construye las matrices. Las opciones se escriben de memoria y deben contrastarse con la ayuda de la versión instalada (`--help`).

```bash
# Cuantificar todas las muestras con el mismo archivo samples.txt que se usó para ensamblar
# (sbatch, un job por muestra o un job array en el clúster)
align_and_estimate_abundance.pl \
  --transcripts trinity_proyecto.Trinity.fasta \
  --seqType fq --samples_file samples.txt \
  --est_method salmon --trinity_mode --prep_reference \
  --thread_count 8 --output_dir quant_salmon

# Variante clásica con RSEM y bowtie2
align_and_estimate_abundance.pl \
  --transcripts trinity_proyecto.Trinity.fasta \
  --seqType fq --samples_file samples.txt \
  --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference \
  --thread_count 8 --output_dir quant_rsem

# Matrices de conteo y de TPM a nivel gen y a nivel transcrito
ls quant_salmon/*/quant.sf > quant_files.txt
abundance_estimates_to_matrix.pl --est_method salmon \
  --gene_trans_map trinity_proyecto.Trinity.fasta.gene_trans_map \
  --quant_files quant_files.txt --name_sample_by_basedir \
  --out_prefix salmon
```

### 3.4 Flujo alternativo, paso a paso con Salmon y R

Es el que se usa en los scripts del curso: se indexa el FASTA, se cuantifica cada muestra y se agrega con `tximport`.

```bash
salmon index -t trinity_proyecto.Trinity.fasta -i salmon_index -k 31 -p 8

# por muestra (job array en sbatch)
salmon quant -i salmon_index -l A \
  -1 muestra_R1.trim.fq.gz -2 muestra_R2.trim.fq.gz \
  -p 8 --seqBias --gcBias -o quants/muestra
```

```r
library(tximport)
m <- read.delim("trinity_proyecto.Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(TXNAME = m[[2]], GENEID = m[[1]])
files <- file.path("quants", muestras, "quant.sf"); names(files) <- muestras
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)   # conteos a nivel gen
```

### 3.5 Nivel transcrito frente a nivel gen

| Nivel | Ventaja | Limitación en de novo |
| --- | --- | --- |
| Gen de Trinity | Más estable ante lecturas ambiguas entre isoformas | El "gen" puede mezclar palólogos o dividir alelos |
| Transcrito (isoforma) | Permite buscar uso diferencial de isoformas | La reconstrucción de isoformas es incompleta, y las lecturas ambiguas dan estimaciones ruidosas |

Para expresión diferencial en un ensamblaje de novo, el estándar práctico es el nivel gen.

### 3.6 Control de calidad de la cuantificación

| Revisión | Qué mirar | Si algo falla |
| --- | --- | --- |
| Tasa de mapeo por muestra | Porcentaje de lecturas mapeadas en `aux_info/meta_info.json` o en el log de Salmon | Una muestra con tasa muy baja respecto de las otras: contaminación, otra especie o librería mal indicada |
| Tipo de librería | `lib_format_counts.json` de Salmon (`-l A` lo detecta) | Si no coincide con lo esperado, revisar el protocolo de librería |
| Correlación entre réplicas | Correlación de log(TPM + 1) entre réplicas de una condición | Réplicas poco correlacionadas: ver PCA y outliers (sesión de modelos estadísticos) |
| Genes de baja expresión | Muchos "genes" con conteos casi nulos | Filtrar antes del análisis (por ejemplo, conteos mínimos en un número mínimo de muestras) |
| Redundancia | Número de "genes" desproporcionado y BUSCO con muchos duplicados | Filtrar por expresión o por isoforma más larga y reindexar |

### 3.7 Errores frecuentes

- **Cuantificar con lecturas sin recortar** cuando el ensamblaje se hizo con lecturas recortadas: se pierde consistencia entre muestras.
- **Usar TPM como entrada de DESeq2 o edgeR:** estos métodos piden conteos (Sesión de modelos estadísticos).
- **Invertir el orden de columnas del mapa gen-transcrito** en `tximport`: Trinity trae gen primero y transcrito después.
- **Cambiar el FASTA después de cuantificar:** los nombres dejan de coincidir; hay que reindexar y volver a cuantificar.
- **Mezclar muestras cuantificadas contra ensamblajes distintos** en una misma matriz.

## Cierre, preguntas y referencias

Trinity ensambla en tres módulos con demandas de recursos distintas; el ensamblaje se juzga con medidas complementarias (ExN50, BUSCO, back-mapping) y se anota con TransDecoder y homología; y se cuantifica cada muestra contra esa referencia única, agregando a "gen" de Trinity.

### Preguntas para discusión (5 min opcionales)

1. ¿Por qué Inchworm y Chrysalis piden memoria y Butterfly pide núcleos, y cómo cambia eso la solicitud en SLURM?
2. Un ensamblaje tiene N50 alto, pero BUSCO muestra muchos duplicados. ¿Qué medida adicional pedirían y qué harían?
3. Solo el 20 % de los transcritos tiene anotación. ¿Es un problema del ensamblaje, de la base de datos o ninguno de los dos?
4. ¿Por qué se cuantifican todas las muestras contra un mismo ensamblaje en lugar de ensamblar cada muestra?
5. Una muestra tiene 45 % de mapeo y el resto 90 %. ¿Qué revisarían primero?

### Referencias

Citadas de memoria o tomadas de las páginas consultadas; conviene verificar volumen y DOI antes de proyectarlas.

- Grabherr MG et al. 2011. Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nature Biotechnology.
- Haas BJ et al. 2013. De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis. Nature Protocols.
- Simão FA et al. 2015. BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics.
- Bryant DM et al. 2017. A tissue-mapped axolotl de novo transcriptome enables identification of limb regeneration factors. Cell Reports (Trinotate).
- Smith-Unna R et al. 2016. TransRate: reference-free quality assessment of de novo transcriptome assemblies. Genome Research.
- Li B, Dewey CN. 2011. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics.
- Patro R et al. 2017. Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods.
- Bray NL et al. 2016. Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology.
- Soneson C, Love MI, Robinson MD. 2015. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research.

### Fuentes consultadas para este documento

- [Grabherr et al., artículo original de Trinity (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3571712) y [copia con la descripción de Inchworm y k = 25](https://scispace.com/pdf/full-length-transcriptome-assembly-from-rna-seq-data-without-55o8f60wxt.pdf).
- [Haas et al., protocolo de Trinity (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3875132).
- [Wiki de Trinity en GitHub](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Home)
- [CUNY HPC](https://wiki.csi.cuny.edu/cunyhpc/index.php/TRINITY), [MetaCentrum](https://wiki.metacentrum.cz/wiki/Trinity), [FASRC Harvard](https://docs.rc.fas.harvard.edu/kb/trinity/) y [St Andrews](https://stab.st-andrews.ac.uk/wiki/index.php/Trinity)
- [NYU Genomics Core](https://learn.gencore.bio.nyu.edu/denovo/transcriptome.html), [OICR](https://bioinformaticsdotca.github.io/rnaseq_2017_tutorial6) y [UFS](https://docs.hpc.ufs.ac.za/training/transcriptomics_tutorial/transcriptomics_tutorial_1/)
- [Eel Pond](https://eel-pond.readthedocs.io/en/latest/4-evaluating-assembly.html)
- [MacManes, preprint](https://www.biorxiv.org/content/10.1101/035642v1.full-text) y [Duarte et al., Bio-protocol](https://bio-protocol.org/en/bpdetail?id=3912)

