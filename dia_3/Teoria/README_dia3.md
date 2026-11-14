## Estructura del Curso

- **Bloque Teórico (1 hora):** Historia del análisis transcriptómico, alineamiento *splice-aware* vs. pseudoalineamiento, métricas de normalización (TPM, CPM, FPKM/RPKM) y preparación de matrices para análisis de Expresión Diferencial (DEGs).
- **Práctica en HPC con Nextflow (3 horas):** Construcción de índices, alineamiento vs. pseudoalineamiento, obtención de matrices de conteo y comparación de estrategias de anotación.

---

## Bloque Teórico (1 Hora)

### 1. Contexto Histórico: Evolución de las Herramientas en RNA-seq

El análisis de transcriptómica ha evolucionado paralelamente al incremento en el rendimiento de la secuenciación de segunda generación (Next-Generation Sequencing, NGS):

- **Primera Generación (2009–2013):** Con la llegada de las plataformas Illumina, surgieron los primeros alineadores con capacidad de empalme como **TopHat** (Trapnell et al., 2009), el cual utilizaba **Bowtie** como motor primario para mapear lecturas no divididas y luego buscaba exones vecinos.
- **Segunda Generación (2013–2016):** La baja velocidad y alta tasa de falsos positivos de TopHat impulsaron el desarrollo de alineadores de segunda generación. **STAR** (Dobin et al., 2013) introdujo el uso de tablas de sufijos (*Suffix Arrays*) para lograr velocidades de alineamiento órdenes de magnitud superiores, mientras que **HISAT** y **HISAT2** (Kim et al., 2015, 2019) optimizaron el uso de memoria RAM mediante la Transformada de Burrows-Wheeler Jerárquica (GFM).
- **Tercera Generación (2016–Presente): Era de la Cuantificación Ligera (Pseudoalineamiento):** La necesidad de procesar miles de muestras en tiempos reducidos derivó en herramientas como **Kallisto** (Bray et al., 2016) y **Salmon** (Patro et al., 2017). Estos métodos omiten el alineamiento base por base a nivel de genoma y cuantifican directamente sobre el transcriptoma utilizando grafos de de Bruijn o $k$-meros.

---

### 2. Alineamiento de Lecturas con Referencia (*Splice-Aware*)

A diferencia del ADN genómico, las lecturas de RNA-seq derivadas de ARNm maduro no contienen intrones. Por ello, los alineadores de transcriptómica deben reconocer sitios de empalme (*splice junctions*) mapeando lecturas a través de exones distantes.

$$\text{Lectura RNA-seq: } \text{[ Exón 1 ]} \text{--------------} \text{[ Exón 2 ]}$$

$$\text{Genoma de Referencia: } \text{======Exón 1======} \dots \text{Intrón} \dots \text{======Exón 2======}$$

#### Fundamentos Algorítmicos

- **STAR (Spliced Transcripts Alignment to a Reference):** Utiliza un algoritmo de búsqueda secuencial de semillas sin comprimir sobre arreglos de sufijos (*Uncompressed Suffix Arrays*). Identifica la coincidencia máxima exacta (MMP - *Maximum Mapped Prefix*) y cose las semillas divididas por intrones.
- **HISAT2 (Hierarchical Indexing for Spliced Alignment of Transcripts):** Basado en la Transformada de Burrows-Wheeler (BWT) y el FM-Index. Emplea dos tipos de índices: un índice global para todo el genoma y miles de pequeños índices locales que cubren regiones genómicas de ~$64,000$ pb, permitiendo una asignación rápida de exones disjuntos con un requerimiento de memoria reducido.

#### Comparativa: HISAT2 vs. STAR vs. Pseudoalineadores (Salmon/Kallisto)

| Característica | HISAT2 | STAR | Salmon / Kallisto |
| :--- | :--- | :--- | :--- |
| **Algoritmo principal** | Hierarchical BWT / GFM-Index | Sequential Seed Search + Suffix Array | Grafos de de Bruijn / Quasi-mapping |
| **Uso de Memoria (RAM)** | **Bajo (~4–8 GB)** | **Alto (~30–32 GB)** | **Muy Bajo (~2–4 GB)** |
| **Velocidad** | Moderada - Alta | **Extremadamente rápido** | **Ultra-rápido (Segundos/Muestra)** |
| **Sensibilidad en *splicing*** | Buena | Excelente (uniones *de novo*) | No aplica (cuantifica transcriptoma) |
| **Output primario** | Archivo BAM genómico | Archivo BAM genómico | Matrices de abundancia (Isoformas/Genes) |
| **Escenario ideal** | Equipos con recursos limitados o datasets masivos | Clusters HPC; detección de variantes/splicing *de novo* | Estimación rápida de expresión diferencial génica/de transcritos |

---

### 3. Unidades de Medida y Normalización en RNA-seq

Para comparar expresión entre muestras, es necesario corregir por dos sesgos técnicos:

1. **Profundidad de secuenciación:** Muestras con mayor número de lecturas tendrán más conteos por gen de forma artificial.
2. **Longitud del transcrito:** Genes más largos acumulan más fragmentos secuenciados que genes cortos expresados al mismo nivel de concentración molecular.

#### Formulación de Métricas

- **RPKM / FPKM (*Reads/Fragments Per Kilobase Million* - Mortazavi et al., 2008):**

$$\text{FPKM} = \frac{\text{Conteo del gen} \times 10^9}{\text{Longitud del gen (pb)} \times \text{Lecturas totales mapeadas}}$$

*Inconveniente histórico:* La suma total de FPKM varía entre muestras debido al orden de normalización, lo que invalida la comparación directa de abundancias relativas (Wagner et al., 2012).

- **TPM (*Transcripts Per Million* - Li et al., 2010):** Normaliza **primero por longitud del gen** y **luego por profundidad de librería**.

$$\text{TPM}_i = \frac{\frac{\text{Conteo}_i}{\text{Longitud}_i}}{\sum_j \left(\frac{\text{Conteo}_j}{\text{Longitud}_j}\right)} \times 10^6$$

*Ventaja:* La suma de TPM en cada muestra **siempre es exactamente** $1,000,000$, lo que permite interpretar TPM como una medida de concentración relativa de ARNm.

- **CPM (*Counts Per Million*):** Normaliza **únicamente por profundidad de librería** (sin corregir por longitud).

$$\text{CPM}_i = \frac{\text{Conteo}_i}{\text{Lecturas totales}} \times 10^6$$

*Uso principal:* Filtrado inicial de genes de baja expresión.

> **Nota:** Las herramientas de expresión diferencial cuantitativa como `DESeq2` y `edgeR` utilizan modelos estadísticos basados en la distribución **Binomial Negativa** y requieren **estrictamente matrices de conteos crudos** (*raw counts*) o sus equivalentes escalados mediante herramientas como `tximport` (Soneson et al., 2015), nunca valores precortados o normalizados como TPM/FPKM.

---

### 4. Preparación de Datos para DEGs en R/Bioconductor

Los paquetes estadísticos dominantes emplean métodos de normalización robustos a la presencia de genes altamente expresados desproporcionadamente:

- **TMM (*Trimmed Mean of M-values* - Robinson & Oshlack, 2010):** Utilizado por `edgeR`. Calcula factores de escalamiento recortando las colas extremas de las razones de expresión (valores M) y de las expresiones absolutas (valores A).
- **Median of Ratios (Anders & Huber, 2010):** Utilizado por `DESeq2`. Calcula una muestra pseudo-referencia mediante el promedio geométrico de cada gen a través de todas las muestras, calculando la mediana de las razones por muestra para estimar el *size factor*.

```r
library(edgeR)
library(DESeq2)

# 1. Carga de matriz de conteos y metadatos
counts <- read.table("gene_counts_matrix.txt", header = TRUE, row.names = 1)
colData <- read.csv("metadata.csv", row.names = 1)

# Validar concordancia de muestras
stopifnot(all(colnames(counts) == rownames(colData)))

# 2. Filtrado de baja expresión (< 1 CPM en la mayoría de las muestras)
keep <- filterByExpr(counts, group = colData$condition)
counts_filtered <- counts[keep, ]

# 3. Normalización por TMM (Trimmed Mean of M-values) en edgeR
dge <- DGEList(counts = counts_filtered, group = colData$condition)
dge <- calcNormFactors(dge, method = "TMM")

# 4. Normalización equivalente por Median of Ratios en DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = colData,
                              design = ~ condition)
dds <- estimateSizeFactors(dds)
```

---

### 📖 Referencias Bibliográficas

1. **Anders, S., & Huber, W.** (2010). Differential expression analysis for sequence count data. *Genome Biology*, 11(10), R106.
2. **Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L.** (2016). Near-optimal probabilistic RNA-seq quantification. *Nature Biotechnology*, 34(5), 525–527.
3. **Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R.** (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1), 15–21.
4. **Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L.** (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature Biotechnology*, 37(8), 907–915.
5. **Li, B., Ruotti, V., Stewart, R. M., Thomson, J. A., & Dewey, C. N.** (2010). RNA-Seq gene expression estimation with or without a genome reference. *BMC Bioinformatics*, 11(1), 400.
6. **Mortazavi, A., Williams, B. A., McCue, K., Schaeffer, L., & Wold, B.** (2008). Mapping and quantifying mammalian transcriptomes by RNA-Seq. *Nature Methods*, 5(7), 621–628.
7. **Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C.** (2017). Salmon provides fast and bias-aware transcript expression estimates. *Nature Methods*, 14(4), 417–419.
8. **Robinson, M. D., & Oshlack, A.** (2010). A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology*, 11(3), R25.
9. **Soneson, C., Love, M. I., & Robinson, M. D.** (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research*, 4, 1521.
10. **Trapnell, C., Pachter, L., & Salzberg, S. L.** (2009). TopHat: discovering splice junctions with RNA-Seq. *Bioinformatics*, 25(9), 1105–1111.
11. **Wagner, G. P., Kin, K., & Lynch, V. J.** (2012). Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples. *Theory in Biosciences*, 131(4), 281–285.
