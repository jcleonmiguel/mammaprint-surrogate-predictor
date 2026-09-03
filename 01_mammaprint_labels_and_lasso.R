############################################################
# Predictor sustitutivo de MammaPrint a partir de RNA-seq de TCGA-BRCA
# PARTE I: etiquetado pronóstico y primer clasificador (LASSO / elastic net)
#
# Autor: José Carlos León Miguel
# Máster Universitario en Bioinformática - Universidad Europea
#
# Objetivo:
# Implementar un pipeline completo de análisis transcriptómico
# que permita:
#   1) seleccionar un subconjunto luminal de muestras BRCA,
#   2) normalizar datos RNA-seq,
#   3) calcular etiquetas pronósticas tipo MammaPrint,
#   4) entrenar un predictor sustitutivo mediante ML.
#
# USO:
#   Rscript 01_mammaprint_labels_and_lasso.R [expresion.tsv] [clinico.tsv] [salida]
#
# Sin argumentos, los datos se leen de data/ y los resultados se escriben en
# output/. Los datos de TCGA no se incluyen en el repositorio (ver README).
# Este script guarda output/results/entorno_parte1.RData, que necesita
# 02_model_comparison.R.
############################################################
# Limpieza del entorno para garantizar reproducibilidad
# y evitar interferencias de objetos previos.
rm(list = ls())
# Fijación de semilla aleatoria para reproducibilidad
# en particiones train/test y procesos de ML.
set.seed(123)
# Rutas de salida: se admite una carpeta alternativa como tercer argumento
# de línea de comandos. Dentro se crean plots/ y results/.
args <- commandArgs(trailingOnly = TRUE)

OUTPUT_DIR  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "output"
PLOTS_DIR   <- file.path(OUTPUT_DIR, "plots")
RESULTS_DIR <- file.path(OUTPUT_DIR, "results")
dir.create(PLOTS_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

plot_path   <- function(f) file.path(PLOTS_DIR, f)
result_path <- function(f) file.path(RESULTS_DIR, f)

# ============================================================
# 0) INSTALACIÓN / CARGA DE PAQUETES
# ============================================================
# Comprobación de la disponibilidad de BiocManager.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Lista de paquetes CRAN necesarios para:
# - manipulación eficiente de datos (data.table)
# - entrenamiento de modelos supervisados (caret)
# - regresión penalizada elastic net (glmnet)
cran_pkgs <- c("data.table", "caret", "glmnet")
# Instalación automática de paquetes CRAN si no están presentes
# en el entorno de ejecución.
for (pkg in cran_pkgs){
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
# Lista de paquetes Bioconductor empleados en el análisis:
# - edgeR: normalización y análisis RNA-seq
# - genefu: firmas transcriptómicas (incluyendo MammaPrint)
# - org.Hs.eg.db / AnnotationDbi: anotación génica
# - mammaPrintData: datasets y objetos auxiliares de la firma
bioc_pkgs <- c("edgeR", "genefu", "org.Hs.eg.db", "AnnotationDbi", "mammaPrintData")
# Instalación automática de paquetes Bioconductor en caso
# de no estar disponibles. Se desactiva la actualización
# global para mantener estabilidad y reproducibilidad.
for (pkg in bioc_pkgs){
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}
# Carga de librerías necesarias para el pipeline.
library(data.table)
library(edgeR)
library(caret)
library(genefu)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(mammaPrintData)

# Cargar datos internos de mammaPrintData 
for (ds in c("buyseRGcy3","buyseRGcy5","glasRGcy3","glasRGcy5")){
  tryCatch(data(list = ds, package = "mammaPrintData"),
           error = function(e) invisible(NULL))
}

# ============================================================
# 1) CARGA DE ARCHIVOS
# ============================================================
# Las rutas se resuelven en este orden: argumento de línea de comandos y,
# si no se pasa, el nombre por defecto dentro de data/.
DATA_DIR <- "data"

ruta_entrada <- function(i, nombre_por_defecto, descripcion) {
  ruta <- if (length(args) >= i && nzchar(args[i])) args[i] else file.path(DATA_DIR, nombre_por_defecto)
  if (!file.exists(ruta)) {
    stop(sprintf(paste0("No se encuentra el archivo %s: '%s'.\n",
                        "  Coloca los datos en '%s/' o pasa la ruta como argumento (ver README)."),
                 descripcion, ruta, DATA_DIR), call. = FALSE)
  }
  normalizePath(ruta)
}

expr_file <- ruta_entrada(1, "BRCA_exp_matrix.tsv",         "de EXPRESIÓN")
clin_file <- ruta_entrada(2, "clinical_info_TCGA-BRCA.tsv", "CLÍNICO")

cat("Matriz de expresión:", expr_file, "\n")
cat("Datos clínicos     :", clin_file, "\n")
# Lectura de datos en formato tabular
expr_df <- fread(expr_file, data.table = FALSE, fill = TRUE, header = TRUE)
# El clínico de TCGA trae una primera columna sin nombre (índice de fila), de modo
# que la cabecera tiene un campo menos que los datos. fread desplazaría los nombres
# de columna; read.delim detecta ese caso y realinea la cabecera correctamente.
clin    <- read.delim(clin_file, header = TRUE, sep = "\t", check.names = FALSE,
                      row.names = NULL, stringsAsFactors = FALSE)
# Comprobación rápida de dimensiones
cat("\nDim expr_df:", dim(expr_df), "\n")
cat("Dim clin:",    dim(clin),    "\n")

# ============================================================
# 2) MATRIZ GENES x MUESTRAS
# ============================================================
# Definir IDs génicos como rownames y eliminar columna identificadora
rownames(expr_df) <- expr_df[[1]]
expr_df <- expr_df[, -1, drop = FALSE]
# Conversión a matriz numérica (genes x muestras)
expr    <- as.matrix(expr_df)
storage.mode(expr) <- "numeric"
# Comprobación básica de estructura
cat("\nDim expr (genes x muestras):", dim(expr), "\n")
cat("Primeros rownames (genes):\n");  print(head(rownames(expr), 6))
cat("Primeros colnames (muestras):\n"); print(head(colnames(expr), 6))

# ============================================================
# 3) DETECTAR COLUMNA SP#### EN CLÍNICO
# ============================================================
# Identificación automática de la columna con IDs de muestras (SP####)
sp_col <- NULL
for (nm in colnames(clin)){
  v <- clin[[nm]]
  if (is.character(v) && any(grepl("^SP\\d+$", v), na.rm = TRUE)){
    sp_col <- nm
    break
  }
}
if (is.null(sp_col)) stop("No encuentro columna con IDs tipo SP#### en el clínico.")
cat("\nColumna ID muestras:", sp_col, "\n")

# ============================================================
# 4) DETECTAR COLUMNA ER Y HER2
# ============================================================
# Función auxiliar para localizar variables binarias clínicas
detect_binary_col <- function(df, target_vals){
  for (nm in colnames(df)){
    v <- df[[nm]]
    if (!is.character(v)) next
    vals <- unique(na.omit(v))
    if (all(target_vals %in% vals) && length(vals) <= 10) return(nm)
  }
  return(NULL)
}
# Detección automática de estado ER y HER2
# Estado ER: se usa la columna estándar de TCGA (er_status_by_ihc) y solo se
# recurre a la detección automática si esa columna no estuviera presente.
er_col <- if ("er_status_by_ihc" %in% colnames(clin)) {
  "er_status_by_ihc"
} else {
  detect_binary_col(clin, c("Positive", "Negative"))
}
her2_col <- if ("her2_status_by_ihc" %in% colnames(clin)) "her2_status_by_ihc" else NULL

if (is.null(er_col)) stop("No encuentro columna ER con valores Positive/Negative.")
cat("Columna ER:", er_col, "| Columna HER2:", ifelse(is.null(her2_col),"no encontrada", her2_col), "\n")
print(table(clin[[er_col]], useNA = "ifany"))

# ============================================================
# 5) FILTRO LUMINAL (ER+ & HER2 != Positive)
# ============================================================
# Definición del subconjunto luminal: ER+ y HER2 no positivo (o NA)
if (!is.null(her2_col)){
  keep_lum <- clin[[er_col]] == "Positive" &
    (is.na(clin[[her2_col]]) | clin[[her2_col]] != "Positive")
} else {
  keep_lum <- clin[[er_col]] == "Positive"
}
# Filtrado clínico y emparejamiento con muestras disponibles en la matriz de expresión
clin_lum  <- clin[keep_lum, , drop = FALSE]
sp_ids    <- unique(na.omit(clin_lum[[sp_col]]))
keep_cols <- colnames(expr) %in% sp_ids
# Control de consistencia: cuántas muestras luminales se retienen y cuántas están en expr
cat("\nMuestras luminales (clínico):", length(sp_ids),
    "| Casadas con expr:", sum(keep_cols), "\n")
# Submatriz de expresión para el subconjunto luminal
expr_lum <- expr[, keep_cols, drop = FALSE]
# Comprobación mínima de tamaño muestral para evitar análisis inestables
if (ncol(expr_lum) < 15) stop("Muy pocas muestras luminales (<15).")
# Exportación de IDs de muestras utilizadas (trazabilidad/reproducibilidad)
writeLines(colnames(expr_lum), result_path("luminal_SP_ids.txt"))

# ============================================================
# 6) NORMALIZACIÓN TMM + logCPM
# ============================================================
# Normalización TMM (edgeR) y transformación a logCPM
dge    <- DGEList(counts = expr_lum)
dge    <- calcNormFactors(dge, method = "TMM")
logCPM <- cpm(dge, log = TRUE, prior.count = 1)
# Verificación de dimensiones tras normalización
cat("\nDim logCPM (genes x muestras):", dim(logCPM), "\n")

# ============================================================
# 6b) CONVERSIÓN ENSEMBL -> HGNC SYMBOL (si fuera necesario)
# ============================================================
# Diagnóstico del tipo de identificador génico (Ensembl vs símbolo)
cat("\nEjemplo rownames(logCPM):\n"); print(head(rownames(logCPM), 8))

if (any(grepl("^ENSG", rownames(logCPM)))){
  cat("Detectados Ensembl IDs -> convirtiendo a HGNC Symbol...\n")
  # Eliminación de sufijos de versión (p.ej. ENSG... .1)  
  ens_clean <- sub("\\..*", "", rownames(logCPM))   # quitar versión .1 .2 etc.
  # Mapeo Ensembl -> SYMBOL 
  sym <- mapIds(org.Hs.eg.db,
                keys      = ens_clean,
                column    = "SYMBOL",
                keytype   = "ENSEMBL",
                multiVals = "first")
  # Filtrado de genes con símbolo válido
  valid <- !is.na(sym) & sym != ""
  cat("Genes con Symbol válido:", sum(valid), "de", length(valid), "\n")
  # Sustitución de IDs por símbolos en mayúsculas 
  logCPM_sym          <- logCPM[valid, , drop = FALSE]
  rownames(logCPM_sym) <- toupper(sym[valid])
  
  # Deduplicar: conservar el gen de mayor varianza por símbolo
  vars_all  <- apply(logCPM_sym, 1, var, na.rm = TRUE)
  sym_names <- rownames(logCPM_sym)
  uniq_syms <- unique(sym_names)
  
  logCPM_dedup <- do.call(rbind, lapply(uniq_syms, function(s){
    rows <- which(sym_names == s)
    if (length(rows) == 1) return(logCPM_sym[rows, , drop = FALSE])
    best <- rows[which.max(vars_all[rows])]
    logCPM_sym[best, , drop = FALSE]
  }))
  rownames(logCPM_dedup) <- uniq_syms
  logCPM <- logCPM_dedup
  
  cat("Dim logCPM tras conversión y deduplicación:", dim(logCPM), "\n")
  cat("Ejemplo rownames tras conversión:\n"); print(head(rownames(logCPM), 8))
  
} else {
  # Aunque ya sean símbolos, asegurar mayúsculas y deduplicar si es necesario
  rownames(logCPM) <- toupper(rownames(logCPM))
  dup <- duplicated(rownames(logCPM))
  if (any(dup)){
    vars_all  <- apply(logCPM, 1, var, na.rm = TRUE)
    sym_names <- rownames(logCPM)
    uniq_syms <- unique(sym_names)
    logCPM_dedup <- do.call(rbind, lapply(uniq_syms, function(s){
      rows <- which(sym_names == s)
      if (length(rows) == 1) return(logCPM[rows, , drop = FALSE])
      best <- rows[which.max(vars_all[rows])]
      logCPM[best, , drop = FALSE]
    }))
    rownames(logCPM_dedup) <- uniq_syms
    logCPM <- logCPM_dedup
    cat("Deduplicados símbolos repetidos. Dim:", dim(logCPM), "\n")
  }
}

# ============================================================
# 7) PCA PLOT
# ============================================================
# PCA sobre muestras (transposición: muestras x genes) con centrado y escalado
X_pca <- t(logCPM)
pca   <- prcomp(X_pca, center = TRUE, scale. = TRUE)
# Exportación del scatter PC1 vs PC2 como análisis exploratorio
png(plot_path("PCA_scatter.png"), width = 1000, height = 700)
plot(pca$x[,1], pca$x[,2], pch = 19, col = "steelblue",
     xlab = "PC1", ylab = "PC2",
     main = "PCA - Subconjunto Luminal (ER+ & HER2 != Pos)")
dev.off()
cat("PCA guardado en plots/PCA_scatter.png\n")

# ============================================================
# 8) MAMMAPRINT (gene70)
# ============================================================
# --- 8a) SYMBOL -> ENTREZ ---
# Mapeo de símbolos HGNC a identificadores Entrez, requeridos por la firma gene70 (genefu)
entrez <- mapIds(org.Hs.eg.db,
                 keys      = rownames(logCPM),
                 column    = "ENTREZID",
                 keytype   = "SYMBOL",
                 multiVals = "first")
# Filtrado a genes con Entrez válido y extracción de matriz correspondiente
keep_e   <- !is.na(entrez) & entrez != ""
X2       <- logCPM[keep_e, , drop = FALSE]
entrez2  <- as.character(entrez[keep_e])

cat("\n[MP] Genes con ENTREZ válido:", nrow(X2), "de", nrow(logCPM), "\n")

# --- 8b) Colapsar duplicados ENTREZ -> matriz muestras x ENTREZ ---
# Si múltiples símbolos mapean al mismo ENTREZ, se colapsan por media para evitar duplicados
u <- unique(entrez2)

# Construcción genes x muestras y posterior transposición a muestras x genes
mat_list <- vector("list", length(u))
for (i in seq_along(u)){
  rows <- which(entrez2 == u[i])
  if (length(rows) == 1){
    mat_list[[i]] <- X2[rows, ]
  } else {
    mat_list[[i]] <- colMeans(X2[rows, , drop = FALSE])
  }
}

X_entrez_gxs           <- do.call(rbind, mat_list)   # genes x muestras
rownames(X_entrez_gxs) <- u
colnames(X_entrez_gxs) <- colnames(X2)

X_entrez_mat           <- t(X_entrez_gxs)            # muestras x genes
cat("[MP] Dim X_entrez_mat (muestras x ENTREZ):", dim(X_entrez_mat), "\n")

# --- 8c) Comprobar solapamiento con panel gene70 ---
# Control de calidad: se exige un solapamiento mínimo con la firma para poder calcular el score
data("sig.gene70", package = "genefu")
panel_genes <- as.character(na.omit(sig.gene70$EntrezGene.ID))
overlap     <- sum(panel_genes %in% colnames(X_entrez_mat))
cat("[MP] Solapamiento con panel MammaPrint:", overlap, "de", length(panel_genes), "genes\n")

if (overlap < 10){
  stop(paste0(
    "Solapamiento insuficiente (", overlap, " genes).\n",
    "Revisa que expr_df tenga Ensembl IDs estándar (ENSG...) o gene symbols HGNC.\n",
    "Primeros rownames: ", paste(head(rownames(expr), 5), collapse=", ")))
}

# --- 8d) Calcular MammaPrint MANUALMENTE ---
# En ausencia de datasets internos de entrenamiento en algunas versiones,
# se implementa una alternativa coherente con el enfoque original: correlación
# de Pearson entre el perfil de cada muestra y el perfil de buen pronóstico.

# Filtrar a genes del panel presentes en la matriz
sig70_sub <- sig.gene70[!is.na(sig.gene70$EntrezGene.ID) &
                          sig.gene70$EntrezGene.ID %in% colnames(X_entrez_mat), ]
sig70_sub$EntrezGene.ID <- as.character(sig70_sub$EntrezGene.ID)

cat("[MP] Genes del panel usados para score:", nrow(sig70_sub), "\n")

# Alineamiento (muestras x genes del panel) y extracción del perfil de buen pronóstico
X_panel <- X_entrez_mat[, sig70_sub$EntrezGene.ID, drop = FALSE]
profile <- sig70_sub$average.good.prognosis.profile

# Score = correlación de Pearson entre cada muestra y el perfil de buen pronóstico
mp_score <- apply(X_panel, 1, function(sample_expr){
  cor(sample_expr, profile, method = "pearson", use = "complete.obs")
})

cat("\n[MP] Resumen scores MammaPrint (correlación):\n")
print(summary(mp_score))

# Etiquetado binario a partir del score (umbral 0)
risk <- ifelse(mp_score >= 0, "LowRisk", "HighRisk")
names(risk) <- rownames(X_entrez_mat)

cat("\n[MP] Clasificación MammaPrint:\n"); print(table(risk, useNA = "ifany"))

# Histograma de scores
png(plot_path("MammaPrint_score_hist.png"), width = 700, height = 500)
hist(mp_score, breaks = 30, col = "steelblue", border = "white",
     main = "Distribución scores MammaPrint",
     xlab = "Correlación con perfil buen pronóstico")
abline(v = 0, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Umbral (0)"), col = "red", lty = 2, lwd = 2)
dev.off()
# Construcción de tabla de etiquetas y exportación para trazabilidad
labels <- data.frame(
  sample_id        = names(risk),
  mammaprint_class = as.character(risk),
  stringsAsFactors = FALSE
)
labels <- labels[!is.na(labels$mammaprint_class), , drop = FALSE]

cat("\n[MP] Muestras clasificadas:", nrow(labels), "\n")
print(table(labels$mammaprint_class))

write.csv(labels, result_path("mammaprint_labels.csv"), row.names = FALSE)

# Visualización: distribución de clases
png(plot_path("MammaPrint_barplot.png"), width = 700, height = 500)
barplot(table(labels$mammaprint_class),
        col  = c("steelblue","tomato"),
        main = "Distribución MammaPrint (gene70)",
        ylab = "Nº muestras", xlab = "Clase de riesgo")
dev.off()

# ============================================================
# 9) MACHINE LEARNING (caret + glmnet)
# ============================================================
# Preparación de matriz de predictores a partir de logCPM (muestras x genes)
X_mp <- t(logCPM)   
# Intersección de muestras con etiquetas disponibles (control de consistencia)
common_ids <- intersect(rownames(X_mp), labels$sample_id)
cat("\n[ML] Muestras comunes logCPM / MammaPrint labels:", length(common_ids), "\n")

if (length(common_ids) < 15)
  stop("Muy pocas muestras comunes para entrenar el modelo (<15).")
# Alineamiento exacto entre predictores y etiquetas
X_mp2      <- X_mp[common_ids, , drop = FALSE]
labels_ok  <- labels[match(common_ids, labels$sample_id), , drop = FALSE]
stopifnot(all(rownames(X_mp2) == labels_ok$sample_id))

# Selección de predictores: top 1000 genes por varianza (reducción dimensional)
vars    <- apply(X_mp2, 2, var, na.rm = TRUE)
top_n   <- min(1000, ncol(X_mp2))
X_ml    <- X_mp2[, order(vars, decreasing = TRUE)[1:top_n], drop = FALSE]

# Dataset final para caret: etiqueta y predictores numéricos
df_ml   <- data.frame(y = factor(labels_ok$mammaprint_class), X_ml,
                      check.names = TRUE)

cat("\n[ML] Clases:\n"); print(table(df_ml$y))

# Verificación de clasificación binaria
n_clases <- nlevels(df_ml$y)
if (n_clases < 2) stop("Solo una clase MammaPrint; no es posible entrenar un clasificador binario.")

# Partición estratificada train/test (70/30)
set.seed(123)
idx      <- createDataPartition(df_ml$y, p = 0.7, list = FALSE)
train_df <- df_ml[ idx, ]
test_df  <- df_ml[-idx, ]
# ============================================================
# 9b) ANÁLISIS DE ESTRUCTURA DEL DATASET
# ============================================================
# Diagnóstico previo al entrenamiento: distribución, NZV, y redundancia por correlación

# ---------- 9.1 Scatter plot: media vs varianza (genes) ----------
gene_means <- colMeans(X_ml, na.rm = TRUE)
gene_vars  <- apply(X_ml, 2, var, na.rm = TRUE)

png(plot_path("MeanVar_scatter.png"), width = 900, height = 700)
plot(gene_means, gene_vars, pch = 19, cex = 0.6,
     xlab = "Media (logCPM) por gen",
     ylab = "Varianza (logCPM) por gen",
     main = "Estructura del dataset: Media vs Varianza (genes)")
dev.off()
cat("MeanVar scatter guardado.\n")

# ---------- 9.2 Boxplot por muestra ----------
png(plot_path("logCPM_boxplot_samples.png"), width = 1200, height = 700)
boxplot(logCPM, outline = FALSE, las = 2, cex.axis = 0.6,
        main = "Distribución logCPM por muestra (TMM + logCPM)",
        ylab = "logCPM")
dev.off()
cat("Boxplot por muestra guardado.\n")

# ---------- 9.3 Density plot (submuestra de 20) ----------
set.seed(123)
sel_samples <- sample(colnames(logCPM), size = min(20, ncol(logCPM)))

png(plot_path("logCPM_density_samples.png"), width = 900, height = 700)
dens_first <- density(logCPM[, sel_samples[1]], na.rm = TRUE)
plot(dens_first,
     main = "Densidades logCPM (submuestra de 20 muestras)",
     xlab = "logCPM", ylab = "Densidad",
     ylim = c(0, max(sapply(sel_samples, function(s)
       max(density(logCPM[, s], na.rm = TRUE)$y))) * 1.1))
for (j in 2:length(sel_samples)) {
  lines(density(logCPM[, sel_samples[j]], na.rm = TRUE),
        col = j)
}
dev.off()
cat("Density plot guardado.\n")

# ---------- 9.4 Variables dummy ----------
# Verificación de que todos los predictores son numéricos (sin codificación dummy)
pred_cols <- setdiff(colnames(train_df), "y")
pred_is_numeric <- sapply(train_df[, pred_cols, drop = FALSE], is.numeric)
n_non_numeric   <- sum(!pred_is_numeric)

dummy_info <- data.frame(
  total_predictors       = length(pred_cols),
  non_numeric_predictors = n_non_numeric
)
write.csv(dummy_info, result_path("dummy_variables_check.csv"), row.names = FALSE)
cat("Variables dummy check guardado. Predictores no numéricos:", n_non_numeric, "\n")

# ---------- 9.5 Varianza cero o casi cero ----------
nzv_metrics <- nearZeroVar(train_df[, pred_cols, drop = FALSE], saveMetrics = TRUE)
n_nzv       <- sum(nzv_metrics$nzv)
write.csv(nzv_metrics, result_path("nearZeroVar_metrics.csv"), row.names = TRUE)
cat("NearZeroVar calculado. Genes con varianza casi cero:", n_nzv, "\n")

# ---------- 9.6 Correlación entre características ----------
# Calculamos sobre los predictores de train_df
cor_mat <- cor(train_df[, pred_cols, drop = FALSE],
               use = "pairwise.complete.obs")

high_cor_idx <- findCorrelation(cor_mat, cutoff = 0.90)
n_high_cor   <- length(high_cor_idx)

cor_summary <- data.frame(
  predictors_total                          = length(pred_cols),
  predictors_highly_correlated_cutoff_0.90  = n_high_cor
)
write.csv(cor_summary, result_path("correlation_summary.csv"), row.names = FALSE)
cat("Correlación calculada. Predictores altamente correlacionados (>0.90):", n_high_cor, "\n")

# Heatmap ilustrativo para un subconjunto de predictores
top50_idx <- 1:min(50, length(pred_cols))
png(plot_path("correlation_heatmap_top50.png"), width = 900, height = 800)
heatmap(cor_mat[top50_idx, top50_idx],
        main   = "Heatmap de correlación (top 50 predictores)",
        symm   = TRUE,
        col    = colorRampPalette(c("blue","white","red"))(100),
        margins = c(8, 8))
dev.off()
cat("Heatmap de correlación guardado.\n")

cat("\n[ESTRUCTURA] Todos los plots y resultados guardados:\n")
cat("  - plots/MeanVar_scatter.png\n")
cat("  - plots/logCPM_boxplot_samples.png\n")
cat("  - plots/logCPM_density_samples.png\n")
cat("  - plots/correlation_heatmap_top50.png\n")
cat("  - results/dummy_variables_check.csv\n")
cat("  - results/nearZeroVar_metrics.csv\n")
cat("  - results/correlation_summary.csv\n")
# Configuración de validación cruzada y métrica ROC en caret
ctrl <- trainControl(
  method          = "cv",
  number          = 5,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# Forzar que caret esté activo (puede ser enmascarado por genefu/mclust)
library(caret)

# Entrenamiento del clasificador penalizado (glmnet) con centrado/escalado
set.seed(123)
model <- caret::train(
  y ~ .,
  data       = train_df,
  method     = "glmnet",
  metric     = "ROC",
  trControl  = ctrl,
  preProcess = c("center","scale"),
  tuneLength = 10
)

cat("\n[ML] Resultado del modelo:\n")
print(model)

# Predicción en test y evaluación por matriz de confusión
pred_class <- predict(model, newdata = test_df)
cm         <- confusionMatrix(pred_class, test_df$y)
cat("\n[ML] Confusion Matrix:\n"); print(cm)

# Guardar resultados
save(model, file = result_path("caret_glmnet_model.RData"))
write.csv(model$results,      result_path("caret_glmnet_tuning_results.csv"), row.names = FALSE)
write.csv(as.data.frame(cm$table), result_path("confusion_matrix.csv"),       row.names = FALSE)

# Curva ROC y AUC sobre el conjunto de test
prob_test  <- predict(model, newdata = test_df, type = "prob")
clase_pos  <- levels(df_ml$y)[1]
roc_df     <- data.frame(
  obs  = test_df$y,
  pred = prob_test[[clase_pos]]
)

if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

png(plot_path("ROC_curve.png"), width = 700, height = 600)
roc_obj <- roc(response  = roc_df$obs,
               predictor = roc_df$pred,
               levels    = rev(levels(df_ml$y)),
               quiet     = TRUE)
plot(roc_obj, col = "steelblue", lwd = 2,
     main = paste0("ROC - glmnet | AUC = ", round(auc(roc_obj), 3)))
dev.off()
cat("\n[ML] AUC en test:", round(auc(roc_obj), 3), "\n")
# ============================================================
# 10) GUARDADO DEL ENTORNO PARA LA PARTE II
# ============================================================
# 02_model_comparison.R reutiliza la misma matriz de predictores, las
# etiquetas y la partición train/test, de modo que los modelos de la
# parte II sean comparables con el LASSO entrenado aquí.
# 'risk' se alinea con las filas de X_ml para que la parte II reciba etiquetas
# y predictores con el mismo número de muestras y en el mismo orden.
risk <- risk[rownames(X_ml)]
stopifnot(length(risk) == nrow(X_ml), !anyNA(risk))

save(X_ml, risk, train_df, test_df, model, cor_mat, high_cor_idx,
     file = result_path("entorno_parte1.RData"))
cat("Entorno guardado en:", result_path("entorno_parte1.RData"), "\n")

# Resumen final de outputs generados
cat("\n============================================================\n")
cat("FIN DEL ANÁLISIS\n")
cat("Resultados guardados en:", normalizePath(OUTPUT_DIR), "\n")
cat("  - results/luminal_SP_ids.txt\n")
cat("  - results/mammaprint_labels.csv\n")
cat("  - results/caret_glmnet_model.RData\n")
cat("  - results/caret_glmnet_tuning_results.csv\n")
cat("  - results/confusion_matrix.csv\n")
cat("  - results/entorno_parte1.RData\n")
cat("  - plots/PCA_scatter.png\n")
cat("  - plots/MammaPrint_barplot.png\n")
cat("  - plots/ROC_curve.png\n")
cat("============================================================\n")
# ============================================================
# SESSION INFO
# ============================================================
sessionInfo()