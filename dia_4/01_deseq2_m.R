# =============================================================================
# 01_deseq2.R - Expresion diferencial con DESeq2 a partir de kallisto
# Uso: correlo por secciones en RStudio, o completo con Rscript 01_deseq2.R
# Necesita: samples.tsv con las columnas sample_id y condition
# =============================================================================

library(tximport)
library(DESeq2)
library(ggplot2)
library(pheatmap)

# -------------------------
# Rutas y parametros
# -------------------------
KALLISTO_DIR <- "kallisto_counts"   # Carpeta con una subcarpeta por muestra
SAMPLES      <- "samples.tsv"
OUT_DIR      <- "resultados_deseq2"
REFERENCIA   <- "control"           # Nivel de 'condition' contra el que se compara
PADJ         <- 0.05
LFC          <- 1

dir.create(OUT_DIR, showWarnings = FALSE)

# -------------------------
# 1. Metadatos
# -------------------------
meta <- read.delim(SAMPLES)
rownames(meta) <- meta$sample_id
meta$condition <- relevel(factor(meta$condition), ref = REFERENCIA)
meta

# -------------------------
# 2. Importar kallisto
# -------------------------
# El transcriptoma del NCBI no trae un mapa transcrito -> gen,
# asi que el analisis se hace a nivel de transcrito (txOut = TRUE).
files <- file.path(KALLISTO_DIR, meta$sample_id, "abundance.tsv")
names(files) <- meta$sample_id
stopifnot(all(file.exists(files)))

txi <- tximport(files, type = "kallisto", txOut = TRUE)
saveRDS(txi, file.path(OUT_DIR, "txi.rds"))

# -------------------------
# 3. DESeq2
# -------------------------
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ condition)

# Filtro: al menos 10 cuentas en tantas muestras como el grupo mas pequeno
keep <- rowSums(counts(dds) >= 10) >= min(table(meta$condition))
dds  <- dds[keep, ]

dds <- DESeq(dds)
saveRDS(dds, file.path(OUT_DIR, "dds.rds"))

# Datos transformados para graficas
vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, file.path(OUT_DIR, "vsd.rds"))

# -------------------------
# 4. Graficas de control de calidad
# -------------------------
# PCA
p_pca <- plotPCA(vsd, intgroup = "condition") + theme_bw()
ggsave(file.path(OUT_DIR, "PCA.pdf"), p_pca, width = 6, height = 5)

# Distancias entre muestras
dist_mat <- as.matrix(dist(t(assay(vsd))))
pheatmap(dist_mat, filename = file.path(OUT_DIR, "distancias_muestras.pdf"),
         width = 6, height = 5)

# -------------------------
# 5. Resultados
# -------------------------
resultsNames(dds)
coef <- resultsNames(dds)[2]          # p. ej. "condition_calor_vs_control"

res <- lfcShrink(dds, coef = coef, type = "apeglm")
saveRDS(res, file.path(OUT_DIR, "res.rds"))

res_df <- as.data.frame(res)
res_df$transcrito <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
write.table(res_df, file.path(OUT_DIR, "DE_todos.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sig <- subset(res_df, padj < PADJ & abs(log2FoldChange) >= LFC)
write.table(sig, file.path(OUT_DIR, "DE_significativos.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
nrow(sig)

# -------------------------
# 6. Figuras de resultados
# -------------------------
# MA plot
pdf(file.path(OUT_DIR, "MA_plot.pdf"), width = 6, height = 5)
plotMA(res, ylim = c(-5, 5))
dev.off()

# Volcano
res_df$grupo <- "NS"
res_df$grupo[which(res_df$padj < PADJ & res_df$log2FoldChange >=  LFC)] <- "Up"
res_df$grupo[which(res_df$padj < PADJ & res_df$log2FoldChange <= -LFC)] <- "Down"

p_vol <- ggplot(na.omit(res_df), aes(log2FoldChange, -log10(padj), color = grupo)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "grey70")) +
  geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed") +
  geom_hline(yintercept = -log10(PADJ), linetype = "dashed") +
  theme_bw()
ggsave(file.path(OUT_DIR, "volcano.pdf"), p_vol, width = 6, height = 5)

# Heatmap de los 50 transcritos mas significativos
top <- head(rownames(sig), 50)
pheatmap(assay(vsd)[top, ], scale = "row", show_rownames = FALSE,
         annotation_col = meta["condition"],
         filename = file.path(OUT_DIR, "heatmap_top50.pdf"), width = 6, height = 7)

message("Listo. Resultados en ", OUT_DIR, "/")