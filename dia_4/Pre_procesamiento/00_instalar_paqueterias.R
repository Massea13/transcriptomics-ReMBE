# =============================================================================
# 00_instalar_paquetes.R - Instala los paquetes de R que se usan EN EL SERVIDOR
# Uso: sbatch instalar_paquetes.sh   (una sola vez)
# En el servidor solo se corre 00_exportar_counts.R, que necesita tximport.
# DESeq2, edgeR, ggplot2 y pheatmap se usan en RStudio, no aqui.
# =============================================================================

.libPaths()   # La primera ruta debe ser ~/R/4.4

paquetes <- c("tximport")

# BiocManager es el instalador de los paquetes de Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

faltan <- paquetes[!sapply(paquetes, requireNamespace, quietly = TRUE)]
if (length(faltan) > 0) {
  BiocManager::install(faltan, update = FALSE, ask = FALSE)
}

# Reporte final
estado <- sapply(paquetes, requireNamespace, quietly = TRUE)
print(data.frame(paquete = paquetes, instalado = estado, row.names = NULL))

if (all(estado)) {
  message("Todo listo.")
} else {
  message("Faltan: ", paste(paquetes[!estado], collapse = ", "))
}