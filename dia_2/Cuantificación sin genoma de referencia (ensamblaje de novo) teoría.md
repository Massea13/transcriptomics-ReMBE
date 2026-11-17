# 2.-Cuantificación sin Genoma de referencia (ensamblaje de novo): teoría

Esta sección enseña a cuantificar la expresión génica cuando no hay genoma ni transcriptoma de referencia: primero se ensambla el transcriptoma con Trinity y después se pseudoalinean las lecturas contra ese ensamblaje.

## Estructura de la sesión (45 min)

| Bloque | Tema | Tiempo | Idea central |
| --- | --- | --- | --- |
| 4| Consideraciones biológicas | 15 min | De novo es una decisión biológica; sus límites (isoformas, quimeras, redundancia) se miden y se mitigan |
| 5 | Pseudoalineamiento: Salmon y Kallisto | 15 min | Cuantificar sin alinear base a base: compatibilidad de k-mers y EM |
| 6 | Estructura de índices | 15 min | Qué archivos hacen falta, cómo se construye el índice y qué errores evitar |

## Flujo general

```text
FASTQ crudo -> QC/trimming (FastQC, fastp)
            -> ensamblaje de novo (Trinity) -> Trinity.fasta + gene_trans_map
            -> QC del ensamblaje (BUSCO, TrinityStats, back-mapping)
            -> indice (salmon index / kallisto index)
            -> cuantificacion por muestra (salmon quant / kallisto quant)
            -> matriz de conteos a nivel gen (tximport)
            -> expresion diferencial (DESeq2 / edgeR)
```

## Tema 1. Consideraciones biológicas (15 min)

El enfoque de novo es necesario cuando las lecturas no pueden anclarse con confianza a un genoma o a un transcriptoma anotado; en cualquier otro caso, una referencia de buena calidad suele dar resultados más limpios.

**Reparto sugerido:** 5 min cuándo usarlo, 7 min limitaciones, 3 min cómo evaluarlas.

### 1.1 ¿Cuándo es necesario el enfoque de novo?

| Situación | Por qué lleva a de novo |
| --- | --- |
| Organismo no modelo sin genoma | No hay contra qué alinear |
| Genoma disponible pero fragmentado o sin anotar | Los exones y empalmes se pierden o se anotan mal |
| Población o cepa divergente de la referencia | Las lecturas con muchos mismatches no alinean y se subestima la expresión |
| Interés en transcritos ausentes de la anotación | Genes nuevos, isoformas, transcritos de origen viral, simbionte o parasito |
| Mezclas de organismos (metatranscriptómica, hospedero más microbiota) | No existe una referencia única |
| Comparación entre especies sin genoma común | Un ensamblaje por especie, o uno conjunto, es la única ruta |

Criterios prácticos para decidir:

- **Calidad de la referencia:** completitud (BUSCO) y calidad de la anotación, no solo que exista.
- **Divergencia:** si un alineamiento de prueba da una tasa baja de lecturas mapeadas, la referencia no sirve para esas muestras.
- **Pregunta biológica:** si interesa lo no anotado, conviene ensamblar aunque haya genoma.
- **Alternativa intermedia:** Trinity guiado por genoma (`--genome_guided_bam`) cuando hay genoma pero la anotación es pobre.
- **Costo:** Trinity es intensivo en memoria; la guía histórica de su documentación es de aproximadamente 1 GB de RAM por millón de lecturas.

### 1.2 Limitaciones del ensamblaje de novo

**Isoformas.** Trinity nombra los contigs como `TRINITY_DN1000_c0_g1_i2` (DN = componente, `g` = "gen", `i` = isoforma). Ese "gen" es un agrupamiento por k-mers compartidos, no un gen anotado.

- Puede unir palólogos o genes parecidos en un mismo "gen".
- Puede separar alelos de un mismo gen en "genes" distintos.
- Con lecturas cortas, la conexión entre exones lejanos es ambigua y las isoformas se reconstruyen de forma incompleta.
- Por eso la cuantificación más fiable es a nivel de "gen" de Trinity, no de isoforma.

**Quimeras.** Son contigs que unen secuencias de transcritos distintos.

- Causas: repeticiones y regiones de baja complejidad, UTR que se solapan entre genes vecinos, y artefactos de la librería (cambio de plantilla durante la PCR o la transcripción inversa).
- Efecto: un solo contig acumula lecturas de dos genes, y la expresión diferencial de uno se contamina con la del otro.
- Cómo sospecharlas: cobertura muy desigual a lo largo del contig, o dos regiones con hits a proteínas distintas en BLAST.

**Redundancia.** El ensamblaje contiene varios contigs casi idénticos para un mismo transcrito.

- Causas: heterocigosis y alelos, variantes de empalme, errores de secuenciación y fragmentación de un transcrito en varios contigs.
- Efecto sobre la cuantificación: las lecturas multi-mapeadas se reparten entre contigs redundantes, lo que diluye los conteos y aumenta la varianza.
- Efecto sobre la interpretación: el número de "genes" se infla (cientos de miles) y no corresponde al número real de genes.

**Otras limitaciones que conviene nombrar:**

- Contaminación (bacterias, hongos, ADN humano) se ensambla igual que el organismo de interés.
- Transcritos poco expresados quedan fragmentados o ausentes; los muy expresados pueden generar muchas variantes de error.
- La calidad depende de la profundidad, la longitud de lectura y de si la librería es estrandada.
- Todas las muestras deben cuantificarse contra el mismo ensamblaje para que sean comparables.

### 1.3 Cómo evaluar y mitigar

| Problema | Indicador | Mitigación habitual |
| --- | --- | --- |
| Fragmentación y completitud | BUSCO (porcentaje de completos), ExN50 | Más profundidad, lecturas más largas |
| Redundancia | BUSCO con porcentaje alto de duplicados, número de "genes" desproporcionado | Isoforma más larga por gen, filtrar por expresión (TPM), agrupar con CD-HIT-EST o Corset |
| Quimeras | Cobertura anómala, hits a proteínas distintas | Revisar a mano los genes clave en IGV, filtrar y re-ensamblar |
| Representación de las lecturas | Tasa de back-mapping de las lecturas al ensamblaje | Esperar valores altos (en la práctica, cerca o por encima de 80 %); valores bajos indican ensamblaje incompleto o contaminación |
| Contaminación | Hits taxonómicos fuera del grupo esperado (BLAST, Kraken2) | Filtrar contigs antes de indexar |

Mensaje para el aula: el N50 por sí solo engaña en transcriptomas, porque un ensamblaje redundante o quimérico puede tener un N50 alto. Se reportan BUSCO, ExN50 y la tasa de back-mapping juntos.

### 1.4 Software de este bloque

```bash
# Trinity (ensamblaje de novo; requiere sbatch en clústar)
Trinity --seqType fq --samples_file samples.txt --CPU 16 --max_memory 200G --output trinity_proyecto

# Estadísticas del ensamblaje
TrinityStats.pl trinity_proyecto.Trinity.fasta > TrinityStats.txt

# Completitud
busco -i trinity_proyecto.Trinity.fasta -l insecta_odb10 -m transcriptome -c 16 -o busco_asm

# Isoforma más larga por gen (reduce redundancia para anotación y QC)
get_longest_isoform_seq_per_trinity_gene.pl trinity_proyecto.Trinity.fasta > longest.fasta

# Redundancia por similitud (usar con criterio: agrupa genes parecidos)
cd-hit-est -i trinity_proyecto.Trinity.fasta -o asm_cdhit95.fasta -c 0.95 -T 16
```

El linaje de BUSCO (`insecta_odb10`) es un ejemplo; hay que elegir el del grupo taxonómico de cada proyecto.

## Tema 2. Pseudoalineamiento: Salmon y Kallisto (15 min)

El pseudoalineamiento estima de qué transcritos puede provenir cada lectura sin alinearla base a base, y por eso cuantifica una muestra completa en minutos y con poca memoria.

**Reparto sugerido:** 5 min la idea, 5 min Salmon y Kallisto lado a lado, 5 min comandos y salidas.

### 2.1 La idea

Un alineador clásico (HISAT2, STAR, bowtie2) busca la posición exacta de cada lectura. Para cuantificar eso no hace falta: basta saber con qué transcritos es compatible la lectura.

1. Se descompone la lectura en k-mers.
2. Cada k-mer se busca en el índice, que devuelve los transcritos que lo contienen.
3. La intersección de esos conjuntos es el **conjunto de compatibilidad** de la lectura.
4. Las lecturas con el mismo conjunto se agrupan en una **clase de equivalencia**, y solo se guarda el conteo por clase.
5. Un algoritmo **EM** reparte las lecturas ambiguas entre los transcritos según su abundancia estimada.

El EM alterna dos pasos hasta converger: asigna cada clase de equivalencia a sus transcritos en proporción a la abundancia actual (paso E) y recalcula las abundancias con esas asignaciones (paso M). Es la razón por la que las lecturas multi-mapeadas se tratan de forma probabilística y no se descartan.

**Por qué encaja con un ensamblaje de novo:** el ensamblaje es solo un FASTA de transcritos, que es justo lo que piden estas herramientas; no se necesita genoma, GTF ni BAM.

La unidad de abundancia que reportan es el TPM, que usa la longitud efectiva del transcrito:

```latex
\mathrm{TPM}_i = \frac{x_i / \tilde{l}_i}{\sum_j x_j / \tilde{l}_j} \times 10^6
```

donde x\_i son las lecturas asignadas al transcrito i y l̃\_i su longitud efectiva (longitud menos la longitud media del fragmento, más uno).

### 2.2 Salmon y Kallisto

| Aspecto | Salmon | Kallisto |
| --- | --- | --- |
| Idea de mapeo | Mapeo con k-mers y alineamiento selectivo de las lecturas candidatas | Pseudoalineamiento con un grafo de De Bruijn del transcriptoma |
| Modos | Mapeo contra índice, o cuantificación desde un BAM alineado al transcriptoma | Solo pseudoalineamiento contra índice |
| Detección del tipo de librería | Sí (`-l A`) | No; se indica con `--rf-stranded` o `--fr-stranded` |
| Corrección de sesgos | Secuencia (`--seqBias`), contenido GC (`--gcBias`) y posición (`--posBias`) | Corrección de sesgo de secuencia (`--bias`) |
| Incertidumbre | Bootstraps o muestras de Gibbs | Bootstraps |
| Salida principal | `quant.sf` | `abundance.tsv` y `abundance.h5` |
| Lectura simple (single-end) | `-r archivo.fq` | Requiere `--single -l <largo> -s <sd>` |
| Uso con ensamblajes de novo | Muy extendido; es el método por defecto de Trinity | Igualmente válido y disponible en el mismo script de Trinity |

En la práctica, ambos dan resultados muy parecidos a nivel de gen. La elección suele depender de las funciones que se necesiten (modelos de sesgo, detección automática del tipo de librería) y de qué herramienta esté instalada en el clúster.

**Detalle importante sobre versiones de Salmon.** El alineamiento selectivo (antes `--validateMappings`) es la estrategia de mapeo por defecto desde la versión 1.0.0, y en las versiones recientes esa bandera aparece como obsoleta y sin efecto ([documentación de Salmon 1.3.0](https://gensoft.pasteur.fr/docs/salmon/1.3.0/salmon.html), [opciones de una versión reciente](https://gannet.fish.washington.edu/Atumefaciens/20200429_cbai_salmon_2020GW_transcript_abundances/program_options.log)). Si se usa una versión anterior a la 1.0, hay que añadir la bandera; en las actuales es innecesaria pero inofensiva.

### 2.3 Salidas y control de calidad

```bash
# --- Salmon ---
salmon quant -i salmon_index -l A \
  -1 muestra_R1.trim.fq.gz -2 muestra_R2.trim.fq.gz \
  -p 8 --seqBias --gcBias --numBootstraps 30 \
  -o quants/muestra

# --- Kallisto ---
kallisto quant -i kallisto.idx -o quants_k/muestra -t 8 -b 30 \
  muestra_R1.trim.fq.gz muestra_R2.trim.fq.gz
# estrandado: agregar --rf-stranded (dUTP, el más común) o --fr-stranded
```

| Herramienta | Archivos de salida | Qué revisar |
| --- | --- | --- |
| Salmon | `quant.sf` (Name, Length, EffectiveLength, TPM, NumReads), `aux_info/meta_info.json`, `logs/salmon_quant.log`, `lib_format_counts.json` | Porcentaje de lecturas mapeadas y el tipo de librería inferido |
| Kallisto | `abundance.tsv` (target\_id, length, eff\_length, est\_counts, tpm), `abundance.h5`, `run_info.json` | Lecturas pseudoalineadas del total (en `run_info.json`) |

Lectura de la tasa de mapeo: si es baja, la causa más probable es un ensamblaje incompleto, contaminación en las lecturas o un tipo de librería mal indicado, no un fallo del cuantificador.

### 2.4 Del transcrito al gen

Con ensamblajes de novo la cuantificación por isoforma es poco fiable (Tema 1), así que se agregan los conteos a "gen" de Trinity con `tximport` y el archivo `gene_trans_map`:

```r
library(tximport)
m <- read.delim("trinity_proyecto.Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(TXNAME = m[[2]], GENEID = m[[1]])   # el mapa de Trinity trae gen primero, transcrito después
files <- file.path("quants", muestras, "quant.sf")
names(files) <- muestras
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)   # "kallisto" para abundance.h5
```

Los objetos `txi` se pasan directamente a `DESeqDataSetFromTximport()`.

### 2.5 Alternativa con los scripts de Trinity

```bash
align_and_estimate_abundance.pl --transcripts trinity_proyecto.Trinity.fasta \
  --seqType fq --samples_file samples.txt --est_method salmon \
  --trinity_mode --prep_reference --output_dir salmon_out

abundance_estimates_to_matrix.pl --est_method salmon \
  --gene_trans_map trinity_proyecto.Trinity.fasta.gene_trans_map \
  --name_sample_by_basedir salmon_out/*/quant.sf
```

## Tema 3. Estructura de índices (15 min)

El índice es una estructura que permite buscar rápido en qué transcritos aparece cada k-mer, y se construye una sola vez a partir del FASTA del ensamblaje; en un proyecto de novo ese FASTA es el único "genoma" de referencia, así que lo que entra al índice define todo lo que se puede cuantificar.

**Reparto sugerido:** 5 min archivos de entrada, 5 min cómo es y cómo se construye el índice, 5 min errores frecuentes y buenas prácticas.

### 3.1 Qué archivos se necesitan

| Archivo | Para qué sirve | De dónde sale |
| --- | --- | --- |
| FASTA del ensamblaje (`*.Trinity.fasta`) | Es el transcriptoma de referencia que se indexa | Salida de Trinity |
| Mapa gen-transcrito (`*.gene_trans_map`) | Agrega la expresión de transcritos a "genes" de Trinity | Lo genera Trinity; si falta, se crea con `get_Trinity_gene_to_trans_map.pl` |
| Lecturas recortadas por muestra (FASTQ) | Se cuantifican contra el índice | Salida de fastp o similar |
| Tabla de muestras | Une cada muestra con su condición, lote y archivos | La prepara cada proyecto |
| Índice (directorio de Salmon o archivo `.idx` de Kallisto) | Estructura de búsqueda de k-mers | `salmon index` o `kallisto index` |

No se necesitan genoma, anotación GTF ni archivos BAM.

El encabezado de un contig de Trinity tiene este aspecto:

```text
>TRINITY_DN1000_c0_g1_i1 len=247 path=[0:0-246]
```

Salmon y Kallisto toman como nombre del transcrito lo que va antes del primer espacio, que aquí es `TRINITY_DN1000_c0_g1_i1`. Ese nombre debe coincidir con la segunda columna del mapa gen-transcrito.

### 3.2 Cómo es un índice por dentro

El índice guarda, para cada k-mer del transcriptoma, en qué transcritos y posiciones aparece. Las dos herramientas lo resuelven con estructuras distintas.

| Aspecto | Salmon | Kallisto |
| --- | --- | --- |
| Estructura base | Grafo de De Bruijn compactado y coloreado (pufferfish), con función hash perfecta mínima | Grafo de De Bruijn del transcriptoma con tabla hash de k-mers |
| Formato en disco | Directorio con varios archivos binarios | Un solo archivo (`.idx`) |
| Tamaño de k | Por defecto 31 (impar) | Por defecto 31 (impar, máximo 31) |
| Duplicados | Las secuencias idénticas se colapsan y quedan registradas en `duplicate_clusters.tsv`; se evita con `--keepDuplicates` | Un k-mer compartido por varios transcritos apunta a todos ellos |
| Dependencia de versión | Un índice de una versión mayor distinta puede no cargarse | Igual: se reconstruye al cambiar de versión |

Contenido típico del directorio de índice de Salmon 1.x (puede variar entre versiones):

```text
salmon_index/
  seq.bin, pos.bin, ctable.bin, ctg_offsets.bin   # secuencia y posiciones de los contigs del grafo
  mphf.bin, rank.bin                              # función hash perfecta mínima
  refseq.bin, reflengths.bin, refAccumLengths.bin, complete_ref_lens.bin  # transcritos y longitudes
  duplicate_clusters.tsv                          # transcritos con secuencia idéntica
  info.json, versionInfo.json                     # k usado y versión de Salmon
  pre_indexing.log, ref_indexing.log              # registros de construcción
```

No se editan ni se copian por partes: el directorio completo es el índice.

**Elección de k.** Un k menor aumenta la sensibilidad con lecturas cortas o con divergencia entre lecturas y transcritos; un k mayor aumenta la especificidad. La documentación de Salmon indica que, con alineamiento selectivo, un k menor puede mejorar aún más la sensibilidad ([documentación de Salmon](https://salmon.readthedocs.io/en/latest/salmon.html)). Para lecturas de 100 a 150 pb, el valor por defecto es adecuado.

### 3.3 Cómo se construye

```bash
# Salmon  (requiere sbatch en clúster para ensamblajes grandes)
salmon index -t trinity_proyecto.Trinity.fasta -i salmon_index -k 31 -p 8

# Kallisto
kallisto index -i kallisto.idx -k 31 trinity_proyecto.Trinity.fasta

# Mapa gen-transcrito, si no existe
get_Trinity_gene_to_trans_map.pl trinity_proyecto.Trinity.fasta > trinity_proyecto.Trinity.fasta.gene_trans_map
```

El índice se construye una vez por ensamblaje y se reutiliza para todas las muestras. Si el FASTA cambia (se filtra, se agrupa o se reensambla), hay que reconstruirlo y recuantificar todo.

**Decoys.** En proyectos con genoma, Salmon permite añadir secuencias "señuelo" (genoma) para que las lecturas que no vienen del transcriptoma no se asignen a él por error. En de novo puro no hay genoma, pero si se conoce el genoma de un contaminante o del hospedero, puede usarse como señuelo.

### 3.4 Errores frecuentes y buenas prácticas

- **Indexar un ensamblaje sin depurar.** Contaminantes, ARN ribosomal y contigs diminutos entran al universo de cuantificación; conviene filtrarlos antes.
- **Cambiar el FASTA después de cuantificar.** Los nombres de `quant.sf` dejan de coincidir con el mapa gen-transcrito.
- **Confundir el orden de columnas.** El mapa de Trinity trae gen primero y transcrito después; `tximport` espera lo contrario (transcrito primero).
- **Usar un índice de otra versión.** Anotar en un archivo de texto la versión de Salmon o Kallisto, el valor de k y la ruta del FASTA.
- **Perder la trazabilidad.** Guardar juntos el FASTA, el mapa gen-transcrito, el índice y el registro de comandos usados.
- **No revisar el log de indexado.** Ahí aparecen los duplicados colapsados y las secuencias descartadas.

### 3.5 Resumen en una línea por herramienta

| Paso | Salmon | Kallisto |
| --- | --- | --- |
| Indexar | `salmon index -t asm.fasta -i idx -k 31` | `kallisto index -i asm.idx -k 31 asm.fasta` |
| Cuantificar | `salmon quant -i idx -l A -1 R1 -2 R2 -o out` | `kallisto quant -i asm.idx -o out R1 R2` |
| Resultado | `out/quant.sf` | `out/abundance.tsv` |
| Hacia R | `tximport(..., type = "salmon")` | `tximport(..., type = "kallisto")` |

## Cierre, preguntas y referencias

Se elige de novo cuando la referencia no sirve para la pregunta; sus límites se verifican con BUSCO, ExN50 y la tasa de mapeo; y se cuantifica contra el FASTA del ensamblaje con Salmon o Kallisto, agregando a "gen" de Trinity.

### Preguntas para discusión (5 min opcionales)

1. Un grupo tiene un genoma de una especie hermana. ¿Alinear contra él o ensamblar de novo? ¿Qué prueba rápida harían para decidir?
2. Un ensamblaje tiene N50 alto pero BUSCO con muchos duplicados. ¿Qué pasa con la cuantificación y qué harían?
3. La tasa de mapeo de una muestra es 45 % y la del resto es 90 %. Mencionen tres explicaciones posibles.
4. ¿Por qué agregar a nivel de gen de Trinity es más defendible que reportar isoformas?
5. Cambian el FASTA tras filtrar contaminantes. ¿Qué hay que repetir?

### Referencias

- Grabherr MG et al. 2011. Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nature Biotechnology.
- Haas BJ et al. 2013. De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis. Nature Protocols.
- Bray NL et al. 2016. Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology (Kallisto).
- Patro R et al. 2017. Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods.
- Srivastava A et al. 2020. Alignment and mapping methodology influence transcript abundance estimation. Genome Biology.
- Soneson C, Love MI, Robinson MD. 2015. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research (tximport).
- Simao FA et al. 2015. BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics.
- Davidson NM, Oshlack A. 2014. Corset: enabling differential gene expression analysis for de novo assembled transcriptomes. Genome Biology.

### Fuentes consultadas para este documento

- [Documentación de Salmon, versión 1.3.0](https://gensoft.pasteur.fr/docs/salmon/1.3.0/salmon.html)
- https://salmon.readthedocs.io/en/latest/salmon.html)
- (https://gannet.fish.washington.edu/Atumefaciens/20200429_cbai_salmon_2020GW_transcript_abundances/program_options.log)
- [Guía de Trinity en Galaxy](https://toolshed.g2.bx.psu.edu/repos/iuc/trinity/file/tip/README.rst).


