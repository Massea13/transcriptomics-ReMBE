# =============================================================================
# 00_instalar_paquetes.R - Revisa e instala los paquetes del analisis
# Uso en RStudio: abre el script y presiona "Source" (o Ctrl+Shift+S).
# Solo instala lo que falte; si ya esta todo, no hace nada.
# Si R pregunta si quieres crear una libreria personal, responde "yes".
# =============================================================================

paquetes_cran <- c("ggplot2", "pheatmap")
paquetes_bioc <- c("tximport", "edgeR", "limma", "DESeq2", "sva")

# -------------------------
# 1. Usar paquetes ya compilados (binarios) en Linux tipo RHEL
# -------------------------
# Algunos paquetes (curl, XML, png) necesitan librerias del sistema para
# compilarse. Los binarios de Posit Package Manager ya vienen compilados.
os <- readLines("/etc/os-release")
es_rhel <- any(grepl("rhel|centos|fedora|rocky|almalinux", os, ignore.case = TRUE))
version <- sub('VERSION_ID="?([0-9]+).*', "\\1", grep("^VERSION_ID", os, value = TRUE))

if (es_rhel && version %in% c("8", "9")) {
  options(repos = c(CRAN = paste0("https://packagemanager.posit.co/cran/__linux__/rhel",
                                  version, "/latest")))
  # Posit identifica la version de R por este encabezado
  options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
                                  paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))
  message("Usando binarios de Posit para RHEL ", version)
} else {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  message("Usando CRAN normal (se compilara desde el codigo fuente)")
}

# -------------------------
# 2. Instalar lo que falte
# -------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Primero las dependencias que suelen fallar al compilar
base <- c("png", "curl", "XML", "httr")
faltan_base <- base[!sapply(base, requireNamespace, quietly = TRUE)]
if (length(faltan_base) > 0) install.packages(faltan_base, Ncpus = 4)

faltan_cran <- paquetes_cran[!sapply(paquetes_cran, requireNamespace, quietly = TRUE)]
faltan_bioc <- paquetes_bioc[!sapply(paquetes_bioc, requireNamespace, quietly = TRUE)]

if (length(faltan_cran) > 0) install.packages(faltan_cran, Ncpus = 4)
if (length(faltan_bioc) > 0) BiocManager::install(faltan_bioc, update = FALSE, ask = FALSE, Ncpus = 4)

# -------------------------
# 3. Reporte final
# -------------------------
todos  <- c(base, paquetes_cran, paquetes_bioc)
estado <- sapply(todos, requireNamespace, quietly = TRUE)
data.frame(paquete = todos, instalado = estado, row.names = NULL)