Sesión 3: Construcción de Índices y Cuantificación con transcriptoma de Referencia (4 horas)

**Objetivo: Construir los índices necesarios y ejecutar el flujo de cuantificación de lecturas utilizando un transcriptoma/genoma de referencia.**

Teoría (1 hora):
a) Alineamiento y pseudoalineamiento con referencia. Ventajas y desventajas de cada uno.
b) Unidades de expresión y sesgos técnicos.
c) Preparación y normalización para DEGs.


# **¿Qué es un alineamiento con referencia?**
En el análisis de un **transcriptoma**, el alineamiento con referencia consiste en comparar las lecturas de ARN obtenidas por secuenciación con un **genoma de referencia** o con una colección de transcritos conocidos (**transcriptoma de referencia**).Así se estima de qué gen o transcrito proviene cada lectura y luego se pueden cuantificar los niveles de expresión. (Dobin, 2015) (Mekso, 2022). Una alineación o cuantificación inexacta puede generar falsos positivos o falsos negativos en los análisis posteriores, lo que conlleva conclusiones erróneas.

Existen varios métodos para la alineación y cuantificación, como BWA, Salmon, Kallisto y STAR, desarrollados para abordar los desafíos que plantean los datos de secuenciación de alto rendimiento generados por la secuenciación masiva de ARN. Estos métodos emplean diferentes algoritmos para alinear y cuantificar las lecturas de RNA-seq. Cada uno presenta ventajas y limitaciones, dependiendo del diseño experimental y la calidad de los datos de RNA-seq analizados.
 
Los ARN maduros suelen tener partes que corresponden a exones separados por intrones en el genoma. Por eso, al alinearlas contra el genoma, las herramientas deben permitir alineamientos “con saltos” entre exones y pueden ayudar a detectar empalmes o variantes de transcritos. (Trapnell, 2009) (Dobin, 2015)

Si el organismo no tiene una referencia adecuada, o su secuencia difiere mucho de la referencia disponible, el alineamiento puede perder lecturas o afectar la estimación de expresión. (Benjamin, 2014)(Munger, 2014)

A diferencia del ADN genómico, las lecturas de RNA-seq derivadas de ARNm maduro no contienen intrones. Por ello, los alineadores de transcriptómica deben reconocer sitios de empalme (*splice junctions*) mapeando lecturas a través de exones distantes.

![ADN-ARN](Figuras/ADN-ARN.png)

Figura: Sharma et al., 2025

### 🧬 Pregunta rápida: Estrategias de alineamiento en RNA-seq

Si tu objetivo principal en un análisis de RNA-seq es **descubrir nuevas isoformas de empalme (splicing) o analizar regiones intrónicas y no codificantes**, ¿contra qué referencia debes alinear tus lecturas (*reads*)?

<details>
  <summary><b>A)</b> Contra el transcriptoma de referencia</summary>
  <br>
  ❌ <b>Incorrecto.</b> El transcriptoma de referencia está formado únicamente por secuencias de transcritos ya anotados (principalmente exones maduros), por lo que perderías información de intrones, regiones intergénicas o isoformas no anotadas previamente.
</details>

<details>
  <summary><b>B)</b> Contra el genoma de referencia</summary>
  <br>
  ✅ <b>¡Correcto!</b> El genoma de referencia incluye exones, intrones y regiones no codificantes, lo que permite detectar eventos de <i>splicing</i> alternativo nuevos y transcripción fuera de las anotaciones conocidas (Srivastava, 2020).
</details>

<details>
  <summary><b>C)</b> Es indistinto, ambos contienen exactamente las mismas secuencias</summary>
  <br>
  ❌ <b>Incorrecto.</b> Aunque provienen del mismo organismo, el genoma contiene toda la secuencia genómica (incluyendo intrones y regiones no codificantes), mientras que el transcriptoma solo contiene las secuencias de los transcritos procesados y anotados.
</details>

### 1. Estrategias de Mapeo: Del Genoma al Transcriptoma de Referencia

#### ¿Por qué no usamos STAR ni HISAT2 en esta práctica?
> **💡 Pregunta:** *Si en la mayoría de los artículos de RNA-seq se utilizan STAR o HISAT2, ¿por qué en nuestro flujo de trabajo con oyamel (`Abies`) usamos BWA y Kallisto directamente contra un transcriptoma?*

> **💡 Respuesta:** *Las herramientas como **STAR** y **HISAT2** son alineadores *splice-aware* (conscientes del empalme): necesitan el **genoma completo** (FASTA) y un archivo de coordenadas (GTF/GFF) para "partir" las lecturas a través de los bordes exón-intrón. Sin embargo, las coníferas poseen genomas enormes (~20 Gb, casi 7 veces el tamaño del genoma humano) y altamente repetitivos; alinear contra un genoma de ese tamaño requiere más tiempo y capacidad informática.

Al utilizar un **transcriptoma de referencia** (`GCAT_AB-RNA-1.0.16.fa`), las secuencias provienen de ARNm maduro (exones empalmados sin intrones), por lo que empleamos alineadores lineales rápidos (**BWA**) o pseudoalineadores (**Kallisto**).*


# Alineamiento *splice-aware*

En RNA-seq, una lectura puede abarcar dos exones que están separados por un intrón en el genoma. Un alineador *splice-aware* permite que la lectura se ajuste a dos segmentos del genoma y que el espacio entre ellos se interprete como un intrón, en vez de exigir que toda la lectura coincida de forma continua. (Musich, 2021)(Raplee, 2019)

En términos generales, ambos programas buscan segmentos de la lectura que coincidan con el genoma y evalúan si pueden combinarse en un alineamiento que cruce una unión exón-exón; pueden usar anotaciones de empalme conocidas, aunque también buscan uniones nuevas. (Fenn, 2022)(Dobin, 2015)


## Principios de funcionamiento

* **HISAT2** usa índices FM basados en la transformación de Burrows–Wheeler. El diseño jerárquico de HISAT combina un índice del genoma completo, que ubica anclas de la lectura, con índices locales que ayudan a extender los alineamientos; HISAT2 también puede emplear un índice en forma de grafo que incorpora variantes genómicas. (Kim, 2015) (Kim, 2019)

* **STAR** busca segmentos de lectura —*seeds*— en un índice de sufijos no comprimido, los agrupa y los “ensambla” para obtener alineamientos continuos o divididos por intrones. Su método permite identificar uniones de empalme, incluso algunas no anotadas previamente. (Dobin, 2013)


**El alineamiento *splice-aware* ubica la lectura en el “mapa” del genoma, teniendo en cuenta los intrones; el pseudoalineamiento pregunta “¿con qué transcritos de referencia es compatible esta lectura?”. (Dobin, 2013) (Bray, 2016)**


| Criterio | Alineadores *Splice-Aware* (STAR / HISAT2) | Alineamiento Lineal (BWA) | Pseudoalineamiento (Kallisto) |
| :--- | :--- | :--- | :--- |
| **Referencia requerida** | Genoma completo (FASTA) + Anotación (GTF/GFF) | Transcriptoma (`GCAT_AB-RNA-1.0.16.fa`) | Transcriptoma (`GCAT_AB-RNA-1.0.16.fa`) |
| **¿Resuelve intrones?** | Sí, parte las lecturas en bordes exón-intrón | No (asume secuencias continuas sin intrones) | No (asume secuencias continuas sin intrones) |
| **Consumo de RAM en HPC** | Muy alto (30–100+ GB según el genoma) | Bajo (`8G` en SLURM) | Muy bajo (`8G` en SLURM) |
| **Velocidad de ejecución** | Rápida (si hay suficiente RAM) | Moderada (minutos a horas por muestra) | Extremadamente rápida (segundos a minutos) |
| **Salida principal** | Archivo de alineamiento (`SAM`/`BAM`) | Archivo de alineamiento (`SAM`/`BAM`) | Tabla de abundancias (`abundance.tsv`) |

---

#### Supongamos unos escenarios

**BWA (`bwa mem`)** es como un detective que busca coincidencias **letra por letra (base por base)** en todo el transcriptoma mediante un índice comprimido con la Transformada de Burrows-Wheeler (BWT) y un índice FM (*FM-index*). Cuando entrega su expediente completo (`archivo BAM`) incluye las coordenadas exactas de cada lectura, pero cuando dos isoformas con secuencias idénticas, marca a las lecturas como ambiguas (`MAPQ = 0`).

```bash
# 1. Construir índice BWT (se ejecuta una sola vez sobre la referencia)
bwa index /scratch/taller/datos_oyamel/REFERENCE/GCAT_AB-RNA-1.0.16.fa

# 2. Alinear base por base + convertir SAM a BAM + extraer conteos
bwa mem -t 4 "$REF" "${TRIM_DIR}/${SAMPLE}_1P.fq.gz" "${TRIM_DIR}/${SAMPLE}_2P.fq.gz" | \
    samtools view -@ 4 -bS - | \
    samtools sort -@ 4 -o "${OUT_DIR}/${SAMPLE}.sorted.bam"

samtools index "${OUT_DIR}/${SAMPLE}.sorted.bam"
samtools idxstats "${OUT_DIR}/${SAMPLE}.sorted.bam" > "${OUT_DIR}/${SAMPLE}_counts.txt"

```

**Kallisto (kallisto quant)** es un detective que en lugar de ver letra por letra, opta por utilizar fragmentos para realizar su búsqueda. Puede utilizar palabras cortas de 31 letras (k-mers) para construir un Grafo de Bruijn del Transcriptoma (T-DBG). Después agrupa las lecturas en clases de equivalencia (equivalence classes) según el subconjunto de isoformas compatibles con esos k-mers, y aplica el algoritmo estadístico Expectation-Maximization (EM) para estimar probabilísticamente la abundancia de cada transcrito.

```bash
# 1. Construcción del Grafo de Bruijn de k-mers (.idx)
kallisto index -i "${OUT_DIR}/GCAT_AB-RNA-1.0.16.idx" "$REF"

# 2. Pseudoalineamiento y estimación probabilística de abundancias
kallisto quant \
    -i "${OUT_DIR}/GCAT_AB-RNA-1.0.16.idx" \
    -o "${OUT_DIR}/${SAMPLE}" \
    -t 4 \
    "${TRIM_DIR}/${SAMPLE}_1P.fq.gz" \
    "${TRIM_DIR}/${SAMPLE}_2P.fq.gz"
```

> **💡 Pregunta:** Si el objetivo de un proyecto no es medir niveles de expresión génica, sino detectar mutaciones puntuales (SNPs) en las regiones codificantes de árboles sanos vs. dañados, ¿qué herramienta es indispensable?*

> **💡 Respuesta:** *BWA, ya que genera archivos posicionales .sorted.bam base por base necesarios para el llamado de variantes. Kallisto descarta la posición exacta de los nucleótidos.*


### 2. Unidades de Medida y Normalización Intra-muestra



# **Importancia de la normalización por la longitud del transcrito y la profundidad de secuenciación.**


Para comparar expresión entre muestras, es necesario corregir por dos sesgos técnicos:

1. **Profundidad de secuenciación:** Muestras con mayor número de lecturas tendrán más conteos por gen de forma artificial.
2. **Longitud del transcrito:** Genes más largos acumulan más fragmentos secuenciados que genes cortos expresados al mismo nivel de concentración molecular.


El número de lecturas asignadas a un gen depende, entre otras cosas, de **cuánto se secuenció** y de la **longitud del gen o transcrito**: a igual abundancia, un transcrito más largo puede acumular más lecturas; una muestra con mayor profundidad de secuenciación también suele producir más lecturas por gen. Por eso, los conteos crudos no siempre son comparables directamente entre genes o muestras. (Zhao, 2020) (Robinson, 2010)(p2)

La normalización busca reducir esas diferencias técnicas y, para comparaciones entre muestras, también corregir los cambios en la **composición global** de las lecturas. Por ejemplo, si unos pocos genes son extremadamente abundantes en una muestra, pueden ocupar una mayor proporción de las lecturas y hacer que otros genes parezcan menos expresados. (Robinson, 2010)(p2) (p3)

Para una comparación estándar de **expresión diferencial a nivel de gen**, suele partirse de conteos de lecturas y usar un método basado en modelos de conteo, como **DESeq2** —con factores de tamaño estimados por la mediana de razones— o **edgeR**, que ofrece la normalización TMM. Ambos incorporan factores de normalización al análisis estadístico; no es necesario convertir primero los conteos a TPM para hacer la prueba diferencial. (Love, 2014) (Robinson, 2010)

La corrección por longitud sirve principalmente para estimar y comparar abundancias entre **genes o transcritos de distinto tamaño**. En cambio, al comparar el mismo gen entre muestras, la longitud del gen normalmente es constante; por eso, dividir rutinariamente cada conteo por la longitud antes de una prueba diferencial puede no ser apropiado. TPM y RPKM consideran longitud y profundidad, pero sus valores son relativos a la composición de la muestra y no siempre son directamente comparables entre condiciones. (Zhao, 2020) (Robinson, 2010)

***


# Métricas de normalización (TPM, CPM, FPKM/RPKM)
Estas medidas convierten los conteos de lecturas o fragmentos en valores ajustados por **profundidad de secuenciación** y, en algunos casos, por **longitud**. La diferencia clave es que **CPM no corrige por longitud**, mientras que **RPKM, FPKM y TPM sí lo hacen**. (Zhao, 2021)(Aanes, 2014)


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


### Definiciones y fórmulas

Sea $C\_i$ el conteo asignado al gen o transcrito $i$, $L\_i$ su longitud en bases y $N$ el total de lecturas o fragmentos mapeados de la muestra.

| Unidad   | Nombre                                                | Cálculo conceptual                                                                                           | ¿Corrige longitud? |
| -------- | ----------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ | ------------------ |
| **CPM**  | *Counts per million*                                  | $\displaystyle CPM\_i=\frac{C\_i}{N}\times10^6$                                                              | No                 |
| **RPKM** | *Reads per kilobase per million mapped reads*         | $\displaystyle RPKM\_i=\frac{C\_i}{(L\_i/1000)(N/10^6)}$                                                     | Sí                 |
| **FPKM** | *Fragments per kilobase per million mapped fragments* | $\displaystyle FPKM\_i=\frac{F\_i}{(L\_i/1000)(M/10^6)}$                                                     | Sí                 |
| **TPM**  | *Transcripts per million*                             | Primero se calcula $R\_i=C\_i/(L\_i/1000)$; luego $\displaystyle TPM\_i=\frac{R\_i}{\sum\_j R\_j}\times10^6$ | Sí                 |

En **CPM**, el conteo de cada gen se divide por el tamaño de la biblioteca y se expresa por millón; según el análisis, el conteo puede ser de lecturas o de fragmentos. Como CPM no considera longitud, sirve para ajustar el tamaño de la biblioteca, pero no para comparar directamente la abundancia de genes de longitudes muy distintas. (Aanes, 2014)(p2)

**RPKM** divide los conteos por la longitud del gen o transcrito —en kilobases— y por el número de lecturas mapeadas —en millones—. La corrección por longitud permite comparar, dentro de una muestra, unidades de distinto tamaño; la corrección por millón ajusta la profundidad de secuenciación. (Zhao, 2021)

**FPKM** aplica una idea equivalente a RPKM, pero cuenta **fragmentos** en lugar de lecturas, por lo que se asocia especialmente con secuenciación *paired-end*: los dos extremos secuenciados de un mismo fragmento se consideran una sola unidad de fragmento para cuantificarlo. (Zhao, 2021)(p3) La distinción entre FPKM y RPKM está, por tanto, en qué se cuenta; ambos ajustan por longitud y tamaño de biblioteca. 

**TPM** también corrige por longitud, pero cambia el orden de la normalización: primero divide cada conteo por la longitud y después reescala esos valores para que la suma de TPM de las unidades cuantificadas en cada muestra sea un millón. (Zhao, 2021)(p3) Por esa propiedad, la suma de TPM es igual entre muestras, mientras que la suma de RPKM o FPKM no está fijada del mismo modo. 

### ¿Cómo interpretarlas para genes y transcritos?

A **nivel de transcrito**, cada unidad corresponde a una isoforma y su longitud puede ser distinta de la de otras isoformas del mismo gen; por eso, la cuantificación requiere estimar qué lecturas o fragmentos corresponden a cada transcrito. (Zhao, 2021) A **nivel de gen**, las lecturas de varias isoformas pueden resumirse en un solo valor, pero definir la longitud representativa del gen no siempre es sencillo cuando sus isoformas tienen longitudes diferentes. (Zheng, 2011)

Las longitudes utilizadas por los programas no siempre son simplemente la longitud total anotada: en la cuantificación de transcritos pueden emplearse **longitudes efectivas**, que consideran qué fragmentos podrían haber sido observados dadas las características de la secuenciación. (Lee, 2010) Por tanto, los valores calculados por distintas herramientas pueden diferir si usan métodos o estimaciones de longitud distintos. (Zhao, 2021)

### Ventajas y límites prácticos

* **CPM** es sencillo y útil para expresar conteos ajustados por tamaño de biblioteca, pero no normaliza por longitud; no es una medida directa para comparar genes de tamaños diferentes. (Aanes, 2014)
* **RPKM y FPKM** facilitan la comparación relativa entre genes o transcritos dentro de una muestra al ajustar por longitud y profundidad, pero sus valores pueden depender de la composición de la muestra. (Zhao, 2021)(Zhao, 2020)
* **TPM** facilita la lectura de abundancias relativas dentro de cada muestra porque sus valores suman un millón, pero eso **no garantiza** que un valor de TPM sea directamente comparable entre muestras si cambia mucho la composición de los ARN secuenciados. (Zhao, 2021) (Zhao, 2020)
* Ninguna de estas unidades corrige automáticamente todos los sesgos técnicos, como los relacionados con composición de secuencia, protocolo o cambios específicos de longitud entre muestras. (Zheng, 2011)(Mandelboum, 2019)


Supongamos que evaluamos cuánto se leen dos libros en dos bibliotecas diferentes:
* **Libro A (*El Quijote*):** 1,000 páginas vs. **Libro B (*El Principito*):** 100 páginas (**sesgo por longitud del transcrito**).
* **Biblioteca 1:** 40 millones de lecturas vs. **Biblioteca 2:** 20 millones de lecturas (**sesgo por profundidad de secuenciación**).

Los **conteos crudos** no corrigen ninguno de estos dos factores, por lo que no pueden compararse directamente entre genes distintos ni entre muestras distintas.

| Métrica | ¿Corrige Profundidad? | ¿Corrige Longitud? | ¿Dónde se observa en nuestros outputs? | Uso recomendado |
| :--- | :---: | :---: | :--- | :--- |
| **Conteos Crudos** (*Raw / Est. Counts*) | No | No | Columna 3 de `*_counts.txt` (BWA) y `est_counts` de `abundance.tsv` (Kallisto) | Insumo obligatorio para modelos estadísticos de DEGs (`DESeq2`, `edgeR`). |
| **CPM** (*Counts Per Million*) | Sí | No | Calculado en R (`edgeR::cpm()`) | Filtrar baja expresión; comparar un **mismo transcrito** entre muestras. |
| **RPKM / FPKM** | Sí | Sí | Métrica clásica (Cufflinks/RSEM) | Comparaciones históricas dentro de una misma muestra (la suma total varía entre muestras). |
| **TPM** (*Transcripts Per Million*) | Sí | Sí | Columna `tpm` de `abundance.tsv` (Kallisto) | Visualización exploratoria (*heatmaps*) y comparación de proporciones entre muestras (siempre suma `1,000,000`). |

#### Fórmulas de referencia

* **CPM (*Counts Per Million*):**

$$\text{CPM}_i = \frac{C_i}{\sum_{j} C_j} \times 10^6$$

* **FPKM (*Fragments Per Kilobase Million*):**

$$\text{FPKM}_i = \frac{F_i}{\left(\frac{L_i}{10^3}\right) \times \left(\frac{\sum_{j} F_j}{10^6}\right)}$$

* **TPM (*Transcripts Per Million*):**

$$\text{TPM}_i = \frac{C_i / L_i}{\sum_{j} (C_j / L_j)} \times 10^6$$


### 3. Preparación de Datos para Expresión Diferencial (DEGs) en R/Bioconductor

#### El sesgo de composición de la librería (¿Por qué NO usar TPM para calcular DEGs?)

> **💡 Pregunta:** *En un árbol sano (`DPVR1`) tenemos 4 genes con 250,000 TPM cada uno (suma = 1,000,000). En un árbol dañado por ozono (`DPVR6`), los primeros 3 genes NO cambian su expresión biológica real, pero el Gen 4 (respuesta antioxidante) se sobreexpresa masivamente y acapara el 70% de toda la secuenciación. ¿Qué ocurre con los TPM de los genes 1, 2 y 3?*

> **💡 Pregunta:** *Debido a que los TPM están restringidos a sumar exactamente un millón, **los genes 1, 2 y 3 caerán artificialmente de 250,000 a 100,000 TPM cada uno**. Si se aplicaran pruebas estadísticas sobre TPM, parecería que dichos genes se reprimieron cuando en realidad permanecieron constantes.*

#### Flujo de preparación en R (`edgeR` / `DESeq2`)

1. **Importar conteos y metadatos:**
   * **Desde BWA:** Se unifica la columna 3 de los archivos `*_counts.txt` en una matriz de enteros ($G \text{ transcritos} \times N \text{ muestras}$).
   * **Desde Kallisto:** Se utiliza el paquete `tximport` para importar los archivos `abundance.tsv`, conservando los conteos estimados (`est_counts`) y las longitudes efectivas promedio (`eff_length`) como compensación (*offset*) estadística, **sin redondear valores TPM**.
2. **Filtrado de genes con baja expresión:**
   * Se retienen únicamente los transcritos con al menos **1 CPM** en un mínimo de $k$ muestras, donde $k$ corresponde al tamaño del grupo experimental más pequeño (ej. $k = 3$ para contrastar 3 muestras sanas `DPVR1, 3, 4` vs. 3 dañadas `DPVR6, 8, 9`). Esto reduce el ruido y evita penalizaciones excesivas en la corrección por pruebas múltiples (FDR).
3. **Normalización entre muestras:**
   * **TMM (*Trimmed Mean of M-values* en `edgeR`):** Ordena los genes según su cambio relativo ($\log_2$ fold-change o valor $M$), recorta el 30% superior y el 30% inferior (los genes extremos que responden al tratamiento) y promedia el 40% central para estimar el factor de normalización real de cada librería.
   * **Mediana de Razones (*Median of Ratios* en `DESeq2`):** Construye una pseudo-referencia (media geométrica de cada gen entre todas las muestras), divide cada muestra entre esa referencia y toma la **mediana** como factor de tamaño (*size factor*), evitando que unos pocos genes hiper-expresados distorsionen la escala.



Para un análisis de expresión diferencial a nivel de gen, la matriz debe contener **conteos crudos**: genes en filas, muestras en columnas y números de lecturas o fragmentos asignados a cada gen. No se deben introducir TPM, FPKM, CPM ni conteos previamente normalizados en el análisis de DESeq2. (Love, 2015)(p3)

Los metadatos deben tener una fila por muestra e incluir la condición experimental y, si corresponde, variables como lote, tejido o individuo. Es fundamental que cada fila de metadatos coincida con la columna correspondiente de la matriz de conteos. (Love, 2015)(p10) (Law, 2016)(p4)

### 1. Estructura de los archivos

**`counts.csv`** — primera columna con identificadores de gen y una columna por muestra:

| gene_id    | Muestra_1 | Muestra_2 | Muestra_3 |
| ---------- | --------- | --------- | --------- |
| ENSG000001 | 120       | 98        | 143       |
| ENSG000002 | 0         | 4         | 2         |

**`metadata.csv`** — una fila por muestra y una columna que identifica la condición:

| sample_id | condition | batch |
| --------- | --------- | ----- |
| Muestra_1 | control   | lote1 |
| Muestra_2 | control   | lote2 |
| Muestra_3 | tratado   | lote1 |

### 2. Importar, comprobar y ordenar

El siguiente ejemplo lee ambos CSV, usa los identificadores como nombres de fila y reordena los metadatos para que sigan exactamente el orden de las columnas de conteos. DESeq2 requiere conteos no negativos enteros y una descripción de diseño experimental. (Love, 2015)

```r
library(DESeq2)

# Importar archivos
counts_df <- read.csv(
  "counts.csv",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  "metadata.csv",
  stringsAsFactors = FALSE
)

# Usar el ID de muestra como nombre de fila
rownames(metadata) <- metadata$sample_id

# Comprobaciones básicas de identificadores
stopifnot(!anyDuplicated(rownames(counts_df)))
stopifnot(!anyDuplicated(colnames(counts_df)))
stopifnot(!anyDuplicated(rownames(metadata)))

# Verificar que cada muestra de la matriz esté en los metadatos
stopifnot(all(colnames(counts_df) %in% rownames(metadata)))

# Ordenar metadatos en el mismo orden que las columnas de conteos
metadata <- metadata[colnames(counts_df), , drop = FALSE]

# Confirmar la correspondencia exacta
stopifnot(identical(colnames(counts_df), rownames(metadata)))

# Convertir a matriz y comprobar que sean conteos crudos enteros
counts <- as.matrix(counts_df)
storage.mode(counts) <- "numeric"

stopifnot(all(is.finite(counts)))
stopifnot(all(counts >= 0))
stopifnot(all(counts == round(counts)))

storage.mode(counts) <- "integer"
```

**No redondees TPM/FPKM/CPM para convertirlos en conteos:** esas medidas no son conteos crudos y no cumplen el tipo de entrada esperado por el modelo de DESeq2. (Love, 2015)

### 3. Definir la condición y crear el objeto de DESeq2

Conviene establecer explícitamente el grupo de referencia —por ejemplo, el control— y revisar el diseño antes de ejecutar el análisis. Si hay un lote que debe controlarse y no está confundido con la condición, se puede incluir en la fórmula. (Love, 2015)

```r
metadata$condition <- factor(metadata$condition)
metadata$condition <- relevel(metadata$condition, ref = "control")

# Diseño sencillo:
design_formula <- ~ condition

# Si hay un efecto de lote que se desea modelar:
# metadata$batch <- factor(metadata$batch)
# design_formula <- ~ batch + condition

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = design_formula
)

dds
```

La fórmula `~ condition` compara los niveles de la variable `condition`; una fórmula como `~ batch + condition` incluye el lote como factor del diseño. La estructura del diseño debe reflejar cómo se organizaron realmente las muestras. (Love, 2015)

### 4. Revisiones recomendadas antes de continuar

```r
# Dimensiones: genes × muestras
dim(dds)

# Ver muestras y metadatos asociados
colData(dds)

# Revisar distribución de conteos y tamaño de bibliotecas
summary(colSums(counts(dds)))
summary(rowSums(counts(dds)))

# Revisar cuántos genes tienen cero conteos en todas las muestras
sum(rowSums(counts(dds)) == 0)

# Revisar cuántas muestras hay por condición
table(metadata$condition)
```

Si los conteos se obtuvieron desde archivos de alineamiento, comprueba que los nombres de cromosomas de la anotación coincidan con los de la referencia usada para alinearlos; nombres incompatibles pueden impedir asignar lecturas a genes. (Love, 2015)

Los genes con conteos nulos o muy bajos pueden filtrarse antes del análisis, con un criterio adecuado al diseño y la profundidad; DESeq2 también aplica filtrado independiente en los resultados para reducir pruebas poco informativas. (Love, 2015)(Love, 2014)

**En resumen:** importa los conteos crudos y los metadatos por separado, verifica y alinea explícitamente las muestras, define la condición de referencia y el diseño, y solo entonces crea el objeto `DESeqDataSet`. (Love, 2015)

# Filtrado de genes con baja expresión (ej. < 1 CPM en la mayoría de las muestras).
           Normalización: TMM (Trimmed Mean of M-values) y su implementación en edgeR y DESeq2.

## Filtrado de baja expresión y normalización TMM

Para análisis de expresión diferencial (DEGs), el filtrado y la normalización resuelven problemas distintos: el **filtrado** elimina genes con tan poca información que difícilmente se pueden evaluar con fiabilidad; **TMM** ajusta diferencias de composición entre bibliotecas que no se explican solo por el total de lecturas. (Chen, 2016)

### 1. Filtrar genes de baja expresión

Un umbral como **1 CPM en la mayoría de las muestras** puede servir como orientación, pero no es una regla universal. La cantidad de muestras en que se exige superar el umbral debe corresponder al diseño: por ejemplo, a cuántas réplicas hay en el grupo más pequeño. En un ejemplo de edgeR, se retuvieron genes con más de 0,5 CPM en al menos dos bibliotecas porque había dos réplicas por grupo; los autores indican que el criterio debe conservar genes que tengan expresión suficiente en las réplicas de al menos un grupo. (Chen, 2016)

Un criterio demasiado estricto puede excluir genes que se expresan únicamente en una condición, justo los que podrían ser biológicamente interesantes. En un estudio de pastos, filtrar por 1 CPM eliminó todos los transcritos detectados solo en una de las dos subespecies; ese resultado ilustra el riesgo del umbral en ese conjunto de datos, no que 1 CPM sea siempre inadecuado. (Raithel, 2016b)

En general, el filtrado debe basarse en la cantidad de información de los conteos y no en qué genes parecen cambiar según el contraste, para evitar seleccionar genes usando el resultado que se va a probar. (Chen, 2016)

### 2. ¿Qué hace TMM?

TMM (*trimmed mean of M-values*) estima un factor de escala por muestra. Para ello compara cada biblioteca con una biblioteca de referencia, calcula diferencias de expresión logarítmicas —valores $M$— y niveles de abundancia —valores $A$—, descarta valores extremos y promedia los restantes con ponderaciones. La idea central es que la mayoría de los genes comparables no cambian, o que los cambios al alza y a la baja están razonablemente equilibrados. (Robinson, 2010)(p3) (p7)

En edgeR, los factores TMM **no sustituyen los conteos**: se combinan con el tamaño de biblioteca para calcular el tamaño efectivo, que el modelo usa al analizar expresión diferencial. Por eso, no conviene calcular CPM o dividir conteos y después tratar esos valores transformados como si fueran conteos crudos para el modelo. (Chen, 2016)(p7) (Robinson, 2010)(p6)

### 3. Flujo práctico con edgeR

Ejemplo para comparar dos condiciones. `counts` debe ser una matriz de **conteos crudos**, con genes en filas y muestras en columnas; `meta` debe tener una fila por muestra, en el mismo orden.

```r
library(edgeR)

# Asegurar que el orden de metadatos coincide con las columnas de conteos
meta <- meta[colnames(counts), , drop = FALSE]
stopifnot(identical(rownames(meta), colnames(counts)))

meta$condition <- relevel(factor(meta$condition), ref = "control")

# Crear objeto edgeR
y <- DGEList(counts = counts, group = meta$condition)

# Diseño; añade covariables relevantes si las hay
design <- model.matrix(~ condition, data = meta)

# Filtrado sensible al diseño
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]

# TMM: estima factores de normalización
y <- calcNormFactors(y, method = "TMM")

# Continuar con el modelo de expresión diferencial
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

# Ejemplo: condición tratada frente a control
res <- glmQLFTest(fit, coef = "conditiontreated")
topTags(res)
```

El flujo de edgeR permite retirar genes con conteos insuficientes y usa `calcNormFactors()` para estimar los factores TMM; posteriormente el modelo emplea los tamaños de biblioteca efectivos. (Chen, 2016)(p6) (p7) La estimación de dispersión y el ajuste de un modelo binomial negativo forman parte de la evaluación estadística posterior, no de la normalización TMM en sí. (p10) (p11)

### 4. ¿Y si el análisis se hace con DESeq2?

**TMM es la normalización habitual de edgeR; DESeq2 usa por defecto su propio método, median-of-ratios.** Ambos estiman factores de escala para modelar diferencias entre muestras, pero no son el mismo procedimiento ni sus factores tienen exactamente la misma interpretación. (Love, 2014) (Maza, 2016)

Por tanto, lo más sencillo es usar **TMM dentro del flujo de edgeR** o la normalización predeterminada de **DESeq2** dentro de su propio flujo. DESeq2 permite que el usuario proporcione factores de normalización calculados por otros métodos. (Love, 2014)

Si se necesita explícitamente TMM con DESeq2, se puede calcular el tamaño efectivo de biblioteca en edgeR —tamaño de biblioteca × factor TMM— y pasarlo como factores de tamaño reescalados a media geométrica uno:

```r
library(edgeR)
library(DESeq2)

# counts: conteos crudos; meta: metadatos alineados con las columnas
meta <- meta[colnames(counts), , drop = FALSE]
stopifnot(identical(rownames(meta), colnames(counts)))

meta$condition <- relevel(factor(meta$condition), ref = "control")
design <- model.matrix(~ condition, data = meta)

# Filtrar usando el diseño que se analizará
y <- DGEList(counts = counts)
keep <- filterByExpr(y, design = design)
counts_f <- counts[keep, , drop = FALSE]

# Calcular TMM y tamaños efectivos en edgeR
y_f <- DGEList(counts = counts_f)
y_f <- calcNormFactors(y_f, method = "TMM")
effective_lib <- with(y_f$samples, lib.size * norm.factors)

# DESeq2 requiere factores de escala; se normalizan para que
# su media geométrica sea 1
sf_tmm <- effective_lib / exp(mean(log(effective_lib)))

dds <- DESeqDataSetFromMatrix(
  countData = counts_f,
  colData = meta,
  design = ~ condition
)
sizeFactors(dds) <- sf_tmm
dds <- DESeq(dds)
```

Este uso es una alternativa avanzada, no una razón general para preferir TMM con DESeq2: el método predeterminado de DESeq2 estima factores propios y está integrado en su modelo de conteos. (Love, 2014)(p2) (p14)

### 5. Supuestos y cautelas

TMM corrige sesgos de composición cuando unos pocos genes muy abundantes ocupan una fracción grande de las lecturas. Su premisa es que existe un conjunto suficientemente amplio de genes que sirve como referencia: si la mayoría de los genes cambia en la misma dirección o hay un cambio global real en la cantidad de ARN, la normalización basada solo en conteos puede confundir ese cambio biológico con un efecto de escala. (Robinson, 2010)(p2) (p6)

En experimentos con cambios globales grandes de ARN, TMM y otros métodos de escalado pueden no preservar el cambio absoluto; pueden requerirse controles externos, como *spike-ins*, o una estrategia de normalización basada en una referencia experimental. (Chen, 2016)(p7) (Aanes, 2014)(p2)

**Recomendación práctica:** para un experimento convencional de RNA-seq bulk, aplica un criterio de filtrado relacionado con el tamaño de los grupos y usa `calcNormFactors(y, method="TMM")` para estimar factores TMM antes de analizar los conteos con el modelo de edgeR. Evita imponer “≥1 CPM en la mayoría” sin relacionarlo con el tamaño de los grupos, y revisa si tu pregunta biológica trata de cambios relativos entre genes o de cambios globales en la cantidad total de ARN. (Chen, 2016)(p6) (p7)



---

## Referencias

* Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. *Nature Biotechnology*, 34(5), 525–527. https://doi.org/10.1038/nbt.3519
* Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics*, 25(14), 1754–1760. https://doi.org/10.1093/bioinformatics/btp324
* Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8
* Robinson, M. D., & Oshlack, A. (2010). A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology*, 11(3), R25. https://doi.org/10.1186/gb-2010-11-3-r25
* Soneson, C., Love, M. I., & Robinson, M. D. (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1

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

