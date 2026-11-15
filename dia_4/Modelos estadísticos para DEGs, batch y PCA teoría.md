# Modelos estadísticos para DEGs, batch y PCA: teoría

Esta sección explica cómo se decide estadísticamente qué genes cambian entre condiciones a partir de una matriz de conteos: el modelo (binomial negativa), lo que puede distorsionarlo (outliers y lotes) y la herramienta para mirarlo (PCA).

## Estructura de la sesión (45 min)

| Bloque | Tema | Tiempo | Idea central |
| --- | --- | --- | --- |
| 1 | Modelos estadísticos para DEGs: binomial negativa en DESeq2 y edgeR | 15 min | Los conteos de RNA-seq tienen más varianza que Poisson; la binomial negativa lo modela y las dos herramientas estiman la dispersión compartiendo información entre genes |
| 2 | Outliers y fuentes de variación: efecto de lote | 15 min | Antes de interpretar un DEG hay que saber de dónde viene la variación; el lote se modela en el diseño, no se "borra" |
| 3 | Análisis de componentes principales | 15 min | El PCA es la herramienta de control de calidad que muestra si la muestra se agrupa por condición, por lote o por un outlier |

## Flujo de análisis de esta sesión

```text
quant.sf (Salmon) -> tximport (conteos a nivel gen)
   -> filtrado de genes poco expresados
   -> normalizacion (factores de tamaño en DESeq2, TMM en edgeR)
   -> transformacion (vst / logCPM) -> PCA y distancias entre muestras  [Tema 3]
   -> revisar outliers y lotes, ajustar el diseno                        [Tema 2]
   -> modelo binomial negativo y pruebas de hipotesis                    [Tema 1]
   -> tabla de DEGs (padj, log2FC)
```

El orden de exposición (modelo, variación, PCA) es el de la clase; en un análisis real el PCA se hace primero.

## Tema 1. Modelos estadísticos para DEGs: binomial negativa en DESeq2 y edgeR (15 min)

DESeq2 y edgeR modelan los conteos de cada gen con una distribución binomial negativa porque la variación entre réplicas biológicas es mayor que la que permitiría una Poisson; ambas estiman un parámetro de dispersión por gen y lo estabilizan compartiendo información entre genes.

**Reparto sugerido:** 5 min por qué binomial negativa, 6 min DESeq2 y edgeR paso a paso, 4 min cómo elegir y qué revisar.

### 1.1 Por qué binomial negativa

- Los datos son **conteos de lecturas** (enteros, no negativos), no valores continuos, así que un modelo normal sobre los conteos crudos no es adecuado.
- Una **Poisson** supone que la varianza es igual a la media. Eso se cumple para réplicas técnicas, pero no para réplicas biológicas, que varían más (sobredispersión).
- La **binomial negativa (BN)** añade un término cuadrático a la varianza y se puede ver como una Poisson cuya media varía entre réplicas según una gamma.

```latex
\mathrm{Var}(K_{ij}) = \mu_{ij} + \alpha_i \, \mu_{ij}^2
```

donde K\_ij es el conteo del gen i en la muestra j, μ\_ij su media esperada y α\_i la **dispersión** del gen i. Si α = 0, la BN se reduce a una Poisson. En edgeR el mismo parámetro se llama φ, y su raíz cuadrada es el **coeficiente de variación biológica (BCV)**.

La media se descompone en un factor de normalización por muestra y un componente por gen que depende del diseño:

```latex
\mu_{ij} = s_j \, q_{ij}, \qquad \log_2 q_{ij} = \sum_r x_{jr} \, \beta_{ir}
```

donde s\_j es el factor de tamaño de la muestra j, x\_jr la matriz de diseño (condición, lote) y β\_ir los coeficientes, entre ellos el log2 fold change.

### 1.2 El problema de pocas réplicas y la solución compartida

Con 3 réplicas por grupo, la dispersión de un gen estimada solo con sus datos es muy ruidosa. Las dos herramientas hacen lo mismo en esencia: estiman la dispersión de cada gen, ajustan una tendencia entre dispersión y expresión media, y **encogen** la estimación de cada gen hacia esa tendencia (método bayesiano empírico).

| Aspecto | DESeq2 | edgeR |
| --- | --- | --- |
| Distribución | BN con dispersión α por gen | BN con dispersión φ por gen |
| Normalización | Factores de tamaño por mediana de cocientes (`estimateSizeFactors`) | TMM (`calcNormFactors`) |
| Dispersión | Estimación por gen, tendencia ajustada, encogimiento hacia la tendencia | Dispersión común, tendencia y por gen, con verosimilitud ponderada bayesiana empírica (`estimateDisp`) |
| Prueba habitual | Wald para un coeficiente; LRT para comparar modelos | F de cuasi-verosimilitud (`glmQLFit` y `glmQLFTest`); también `glmLRT` y `exactTest` (dos grupos) |
| Fold change | `lfcShrink` (apeglm, ashr) para estabilizar los log2FC de genes poco informativos | `glmTreat` para probar un umbral de log2FC |
| Outliers | Distancia de Cook; reemplazo de conteos extremos si hay 7 o más réplicas | `robust=TRUE` en `estimateDisp` y `glmQLFit` |
| Filtrado | Independiente (automático en `results`) más filtro previo manual | `filterByExpr` antes de normalizar |
| Corrección de múltiples pruebas | Benjamini-Hochberg (`padj`) | Benjamini-Hochberg (`FDR`) |

Supuestos comunes: la mayoría de los genes no cambia entre condiciones (base de los factores de tamaño y de TMM), las muestras son independientes y los conteos entran como enteros, no como TPM. Para conteos de Salmon, `tximport` aporta las longitudes de transcrito como compensación (offset) en lugar de normalizar por longitud a mano.

Orden de magnitud de la dispersión (guía de edgeR, de memoria): BCV cercano a 0,4 en datos humanos, cercano a 0,1 en organismos modelo genéticamente idénticos. En un ensamblaje de novo con individuos de campo, valores altos son esperables.

### 1.3 Receta en DESeq2

```r
library(DESeq2)
# txi viene de tximport (Tema 2 de la sesión anterior)
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ batch + condition)
keep <- rowSums(counts(dds) >= 10) >= 3          # >= 10 conteos en al menos 3 muestras
dds  <- dds[keep, ]
dds  <- DESeq(dds)                                # tamaños, dispersión, ajuste y prueba de Wald

plotDispEsts(dds)                                 # dispersión: puntos por gen, tendencia y valores finales
res    <- results(dds, contrast = c("condition", "tratamiento", "control"), alpha = 0.05)
resLFC <- lfcShrink(dds, coef = "condition_tratamiento_vs_control", type = "apeglm")
summary(res)
plotMA(resLFC)
```

Si el diseño no incluye lote (una sola fecha de procesamiento), se usa `~ condition`. El nivel de referencia se fija con `relevel`.

### 1.4 Receta en edgeR

```r
library(edgeR)
y <- DGEList(counts = txi$counts, group = meta$condition)
keep <- filterByExpr(y, group = meta$condition)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)                           # TMM

design <- model.matrix(~ batch + condition, data = meta)
y   <- estimateDisp(y, design)
plotBCV(y)                                        # BCV por gen y tendencia
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "conditiontratamiento")
topTags(qlf, n = Inf)$table                       # logFC, logCPM, F, PValue, FDR
```

Este ejemplo pasa `txi$counts` sin compensación por longitud, que es lo más sencillo para enseñar. En un análisis real con tximport conviene seguir la receta de su viñeta, que construye un offset a partir de `txi$length`.

### 1.5 Cómo elegir y qué revisar

- Con diseño sencillo y 3 o más réplicas, **DESeq2 y edgeR dan listas de DEGs muy solapadas**; la elección suele depender de la costumbre del grupo y de las funciones deseadas (`lfcShrink` o `glmTreat`).
- En edgeR se recomienda el flujo de cuasi-verosimilitud (`glmQLFit` y `glmQLFTest`) frente a `glmLRT` para experimentos de RNA-seq en bulk ([discusión en el foro de Bioconductor](https://support.bioconductor.org/p/9136614)).
- **Gráficas de diagnóstico:** dispersión (`plotDispEsts` o `plotBCV`) con puntos alrededor de la tendencia; histograma de p-valores con pico cerca de cero y resto aproximadamente plano; MA plot centrado en cero.
- **Interpretación:** el `padj` controla la proporción esperada de falsos descubrimientos entre los genes declarados; el log2FC mide tamaño de efecto. Se reportan ambos, y el umbral de log2FC se decide antes de mirar los resultados.
- **En de novo:** el "gen" es un agrupamiento de Trinity (Sesión anterior), y un contig quimérico o redundante puede producir DEGs artificiales. Los DEGs más llamativos conviene revisarlos en IGV o contra una base de datos de proteínas.

## Tema 2. Outliers y fuentes de variación: efecto del lote (15 min)

Un DEG solo es interpretable si se sabe qué parte de la variación viene de la biología de interés y qué parte de outliers o de procesamiento; el lote de secuenciación se incorpora al modelo estadístico en lugar de eliminarse de los datos, y solo es corregible si no está confundido con la condición.

**Reparto sugerido:** 4 min fuentes de variación, 4 min outliers, 7 min lote y su corrección.

### 2.1 Fuentes de variación

| Tipo | Ejemplos | ¿Se puede modelar? |
| --- | --- | --- |
| Biológica de interés | Tratamiento, tejido, genotipo | Es lo que se quiere estimar (`condition`) |
| Biológica de fondo | Sexo, edad, individuo, tamaño | Sí, si se registró (covariable o bloque) |
| Técnica conocida | Lote de extracción, fecha de librería, carril o flowcell, operador | Sí, como factor `batch` |
| Calidad de la muestra | Integridad del ARN (RIN), contaminación, porcentaje de ARNr, profundidad | Sí, como covariable, o se excluye la muestra |
| Técnica desconocida | Factores no registrados | Se estiman (variables sustitutas) |
| Específica de de novo | Muestras con tasa de mapeo muy distinta, contigs redundantes o quiméricos | Se revisa en el QC del ensamblaje y de la cuantificación |

### 2.2 Detección de outliers

Un outlier puede ser una **muestra completa** o un **conteo individual** de un gen.

**Muestras.** Se buscan con herramientas que miran todo el transcriptoma a la vez:

1. **PCA** (Tema 3): una muestra aislada, lejos de sus réplicas en PC1 o PC2.
2. **Mapa de calor de distancias** entre muestras (sobre datos transformados con vst) y **agrupamiento jerárquico**.
3. **Métricas de QC:** tasa de mapeo de Salmon, número total de lecturas, RIN, porcentaje de lecturas asignadas.
4. **Correlación entre réplicas** de la misma condición.

**Genes.** DESeq2 calcula la distancia de Cook para marcar conteos extremos. Cuando hay 7 o más réplicas, reemplaza el conteo extremo y reajusta el modelo; el valor por defecto de `minReplicatesForReplace` es 7, y se desactiva con `Inf` ([manual de DESeq2](https://dk-01.installer.hardenedbsd.org/bioconductor/packages/2.14/bioc/manuals/DESeq2/man/DESeq2.pdf), [uso en una consulta del foro de Bioconductor](https://support.bioconductor.org/p/9136684)). Con menos réplicas, el gen se marca con `padj = NA` en lugar de reemplazarse. En edgeR, la opción `robust=TRUE` reduce el peso de genes con comportamiento extremo.

```r
# Distancia de Cook por muestra (DESeq2)
par(mar = c(8, 5, 2, 2))
boxplot(log10(assays(dds)[["cooks"]]), las = 2)

# Distancias entre muestras
vsd <- vst(dds, blind = TRUE)
d   <- dist(t(assay(vsd)))
pheatmap::pheatmap(as.matrix(d), clustering_distance_rows = d, clustering_distance_cols = d)
```

**Qué hacer con un outlier de muestra.** Excluirlo solo con una razón técnica documentada (RIN bajo, contaminación, error de etiquetado, profundidad extrema), no por incomodar al resultado. Se reporta qué muestra se quitó y por qué, y se vuelve a correr todo el análisis. Si la muestra sale rara pero sin causa técnica, se analiza con y sin ella y se informa la diferencia.

### 2.3 Efecto de lote

Un **efecto de lote** es una diferencia sistemática entre grupos de muestras procesadas por separado (distinto día, kit, carril, secuenciador) que no tiene relación con la biología de interés. Se reconoce en el PCA o en el agrupamiento porque las muestras se separan por lote y no por condición.

**Lo primero es el diseño experimental, no la corrección:**

| Diseño | Ejemplo | ¿Corregible? |
| --- | --- | --- |
| Balanceado | Cada lote contiene control y tratamiento en proporción similar | Sí, el lote se incluye en el modelo |
| Parcialmente desbalanceado | Un lote tiene más tratamiento que control | Sí, con menor potencia |
| Confundido | Todo el control en el lote 1 y todo el tratamiento en el lote 2 | No: condición y lote son indistinguibles |

Buenas prácticas de diseño: repartir las condiciones entre lotes y carriles, aleatorizar el orden de procesamiento y registrar todos los metadatos del procesamiento.

**Estrategias de corrección:**

| Estrategia | Cómo | Cuándo usarla | Precauciones |
| --- | --- | --- | --- |
| Incluir el lote en el modelo | `~ batch + condition` en DESeq2; `model.matrix(~ batch + condition)` en edgeR | Es la opción por defecto para pruebas de DEGs | No modifica los datos; necesita réplicas de cada condición en cada lote |
| Quitar el lote solo para visualizar | `limma::removeBatchEffect` sobre datos vst o logCPM | Gráficas (PCA, mapas de calor) | No usar la matriz corregida como entrada del test de DEGs |
| Corregir los conteos | `sva::ComBat_seq` | Cuando se necesita una matriz de conteos corregida para otro análisis | No volver a incluir el lote en el diseño con esos conteos; revisar sobrecorrección |
| Estimar factores ocultos | `sva::svaseq`, `RUVSeq` | Variación no registrada | Las variables sustitutas entran como covariables del diseño; revisar que no absorban la señal biológica |

```r
# 1) Lote en el modelo (la opción recomendada para DEGs)
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ batch + condition)

# 2) Solo para graficar: quitar el lote de los datos transformados
vsd <- vst(dds, blind = TRUE)
mat <- limma::removeBatchEffect(assay(vsd), batch = vsd$batch,
                                design = model.matrix(~ condition, colData(vsd)))

# 3) Conteos corregidos con ComBat-seq
library(sva)
adj <- ComBat_seq(as.matrix(txi$counts), batch = meta$batch, group = meta$condition)

# 4) Variables sustitutas desconocidas
mod  <- model.matrix(~ condition, meta)
mod0 <- model.matrix(~ 1, meta)
norm <- counts(dds, normalized = TRUE)
norm <- norm[rowMeans(norm) > 1, ]
svobj <- svaseq(norm, n.sv = 2, mod, mod0)       # n.sv puede estimarse con num.sv()
# luego: design = ~ SV1 + SV2 + condition, tras añadir svobj$sv a colData
```

Trampas frecuentes:

- **Corregir y después incluir el lote en el diseño con los mismos datos** cuenta el lote dos veces.
- **Quitar el lote de una matriz y pasarla al test** da p-valores demasiado optimistas, porque el modelo ignora los grados de libertad usados en la corrección.
- **Corrección con lote confundido:** la herramienta se ejecuta, pero elimina también el efecto del tratamiento.
- **Un solo lote por condición con una sola réplica:** no hay forma de estimar el efecto del lote.

**Después de corregir,** se repite el PCA: las muestras deben agruparse por condición y no por lote, y el número de DEGs no debe depender drásticamente de la estrategia (si cambia mucho, el resultado es frágil).

## Tema 3. Análisis de componentes principales, PCA (15 min)

El PCA resume en pocas dimensiones la variación entre muestras de miles de genes, y en RNA-seq sirve sobre todo para control de calidad: muestra si las réplicas se agrupan, si el lote o un outlier domina y si la condición separa a las muestras.

**Reparto sugerido:** 5 min la idea, 4 min cómo se calcula bien en RNA-seq, 6 min interpretación y código.

### 3.1 La idea

- Cada muestra es un punto en un espacio con una dimensión por gen (decenas de miles).
- El PCA busca nuevos ejes, los **componentes principales (PC)**, que son combinaciones lineales de genes, ortogonales entre sí, ordenados de mayor a menor varianza explicada.
- **PC1** es la dirección con más varianza entre muestras, **PC2** la siguiente que no esté correlacionada con PC1, y así sucesivamente. Con n muestras hay como máximo n−1 componentes con varianza.
- Cada muestra tiene una **puntuación** (coordenada) en cada PC; cada gen tiene una **carga** (cuánto contribuye a ese PC).
- Se calcula con la descomposición en valores singulares (SVD) de la matriz centrada de genes por muestras, o con los vectores propios de su matriz de covarianza.

```latex
\text{varianza explicada por PC}_k = \frac{\lambda_k}{\sum_m \lambda_m} \times 100\,\%
```

donde λ\_k es el valor propio del componente k.

El PCA es **descriptivo y no supervisado**: no usa las etiquetas de condición y no contrasta hipótesis. Un patrón en el PCA orienta, pero no sustituye a la prueba de DEGs.

### 3.2 Cómo calcularlo bien en RNA-seq

| Decisión | Recomendación | Por qué |
| --- | --- | --- |
| Datos de entrada | Conteos transformados: `vst` o `rlog` en DESeq2; `cpm(y, log = TRUE)` en edgeR | Los conteos crudos tienen varianza que crece con la media; sin transformar, unos pocos genes muy expresados dominan el PCA |
| Genes usados | Los más variables (`plotPCA` usa por defecto 500; se puede cambiar con `ntop`) | Reduce ruido de genes poco informativos ([ejercicio de Cornell sobre `plotPCA`](https://biohpc.cornell.edu/doc/RNA-Seq-2020-exercise2.html)) |
| Centrado y escalado | Centrar por gen; en general no escalar | Escalar da el mismo peso a genes de ruido y a genes con señal |
| `blind` en `vst` | `blind = TRUE` para QC sin influencia del diseño; `blind = FALSE` si se quiere que la transformación use el diseño | Con `blind = FALSE` los grupos pueden verse más separados de lo que están en los datos crudos; declarar cuál se usó |
| Filtrado previo | Quitar genes con conteos muy bajos | Los genes casi sin conteos aportan solo ruido |

```r
library(DESeq2); library(ggplot2)
vsd <- vst(dds, blind = TRUE)

# 1) Atajo de DESeq2 (500 genes más variables)
plotPCA(vsd, intgroup = c("condition", "batch"))

# 2) Versión con control total
mat <- assay(vsd)
top <- order(rowVars(mat), decreasing = TRUE)[1:500]
pca <- prcomp(t(mat[top, ]), center = TRUE, scale. = FALSE)
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)      # % de varianza por PC

df <- data.frame(pca$x[, 1:2], colData(vsd)[, c("condition", "batch")])
ggplot(df, aes(PC1, PC2, color = condition, shape = batch)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", pct[1], "%")) + ylab(paste0("PC2: ", pct[2], "%")) +
  theme_bw()

# 3) Gráfica de sedimentación y genes con mayor carga en PC1
barplot(pct, names.arg = paste0("PC", seq_along(pct)), ylab = "% varianza")
head(sort(abs(pca$rotation[, 1]), decreasing = TRUE), 20)

# edgeR (equivalente): MDS sobre los genes con mayores diferencias
limma::plotMDS(y, col = as.integer(factor(meta$condition)))
```

### 3.3 Cómo leer un PCA

| Patrón observado | Interpretación probable | Qué hacer |
| --- | --- | --- |
| Réplicas juntas y condiciones separadas en PC1 o PC2 | La condición es la fuente principal de variación | Continuar con la prueba de DEGs |
| Muestras separadas por lote y no por condición | Efecto de lote dominante | Incluir `batch` en el diseño (Tema 2); revisar si está confundido |
| Una muestra lejos de todas | Posible outlier técnico | Revisar RIN, tasa de mapeo y conteos; decidir con criterio documentado |
| Condiciones mezcladas y sin estructura | Efecto pequeño, mucha variación biológica o poca potencia | Revisar réplicas y covariables; esperar pocos DEGs |
| PC1 explica un porcentaje muy alto y separa dos grupos | Puede ser una diferencia biológica fuerte o un lote | Colorear por cada metadato disponible antes de concluir |
| Un PC secundario se asocia con RIN, profundidad o fecha | Fuente de variación técnica | Incluirla como covariable |

Buenas prácticas de lectura:

- **Colorear por varias variables** (condición, lote, sexo, RIN) y mirar también PC3 y PC4, no solo PC1 y PC2.
- **Mirar el porcentaje de varianza de cada eje:** dos puntos separados en un PC que explica 3 % pueden ser irrelevantes.
- **Con pocas muestras** (por ejemplo, 6), el PCA es muy sensible a una sola muestra; no conviene sobreinterpretar.
- **El PCA no es una prueba:** que dos grupos se separen no demuestra que haya DEGs, y que no se separen no demuestra que no los haya.
- **En de novo:** una muestra con tasa de mapeo muy distinta puede quedar aparte por razones de ensamblaje y no biológicas.

### 3.4 Actividad corta sugerida (opcional, 3 min)

Mostrar tres PCA en la pantalla (agrupado por condición, agrupado por lote, con un outlier) y preguntar al grupo qué decisión tomarían en cada caso antes de correr la prueba de DEGs.

## Cierre, preguntas y referencias

La binomial negativa modela la sobredispersión de los conteos; el PCA y las distancias entre muestras revelan outliers y lotes antes de la prueba; y el lote se incluye en el diseño siempre que no esté confundido con la condición.

### Preguntas para discusión (5 min opcionales)

1. ¿Por qué una Poisson subestima la variabilidad entre réplicas biológicas y qué consecuencia tiene sobre los falsos positivos?
2. Con 3 réplicas por grupo, ¿por qué conviene encoger la dispersión de cada gen hacia una tendencia?
3. Todo el grupo control se secuenció en un carril y todo el tratamiento en otro. ¿Qué puede y qué no puede concluirse?
4. El PCA separa dos grupos en PC1 (62 %) pero no coinciden con la condición. ¿Qué metadatos revisarían primero?
5. ¿Por qué no se debe usar la matriz corregida con `removeBatchEffect` como entrada del test de DEGs?

### Referencias

- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology.
- Robinson MD, McCarthy DJ, Smyth GK. 2010. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics.
- McCarthy DJ, Chen Y, Smyth GK. 2012. Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research.
- Chen Y, Lun ATL, Smyth GK. 2016. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. F1000Research.
- Robinson MD, Oshlack A. 2010. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology (TMM).
- Zhu A, Ibrahim JG, Love MI. 2019. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics (apeglm).
- Leek JT et al. 2012. The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Bioinformatics.
- Zhang Y, Parmigiani G, Johnson WE. 2020. ComBat-seq: batch effect adjustment for RNA-seq count data. NAR Genomics and Bioinformatics.
- Risso D et al. 2014. Normalization of RNA-seq data using factor analysis of control genes or samples. Nature Biotechnology (RUVSeq).
- Soneson C, Love MI, Robinson MD. 2015. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research (tximport).

### Fuentes consultadas para este documento

- [Manual de DESeq2 (versión antigua)](https://dk-01.installer.hardenedbsd.org/bioconductor/packages/2.14/bioc/manuals/DESeq2/man/DESeq2.pdf)
- [Consulta en el foro de Bioconductor sobre outliers en DESeq2](https://support.bioconductor.org/p/9136684)
- [Ejercicio de RNA-seq de Cornell](https://biohpc.cornell.edu/doc/RNA-Seq-2020-exercise2.html)
- [Discusión en el foro de Bioconductor sobre glmLRT y glmQLFTest](https://support.bioconductor.org/p/9136614)
- [Guía de edgeR en seqtec](https://seqtec.readthedocs.io/en/latest/tools/differential-expression/edger.html)


