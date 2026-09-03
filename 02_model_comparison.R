# =============================================================================
# Predictor sustitutivo de MammaPrint a partir de RNA-seq de TCGA-BRCA
# PARTE II: comparación de modelos de clasificación
#
# Autor: José Carlos León Miguel
# Máster Universitario en Bioinformática - Universidad Europea
# =============================================================================
#
# DESCRIPCIÓN GENERAL DEL SCRIPT:
# Este script entrena varios modelos de clasificación binaria (alto/bajo riesgo
# según MammaPrint) utilizando datos de expresión génica de TCGA-BRCA.
# Los pasos principales son:
#   1. Carga del entorno previo (Entorno de la actividad anterior)
#   2. Inspección y preparación de los datos
#   3. Particionado train/validación
#   4. Entrenamiento con validación cruzada de:
#      - LASSO (glmnet)
#      - Random Forest (rf)
#      - SVM lineal (svmLinear)
#      - Gradient Boosting Machine (gbm)
#   5. Comparación de modelos
#   6. Evaluación final sobre set de validación (matriz de confusión + AUC)
#   7. sessionInfo() para reproducibilidad
#
# USO:
#   Rscript 02_model_comparison.R [entorno_parte1.RData] [carpeta_salida]
#
# Requiere el entorno generado por 01_mammaprint_labels_and_lasso.R. Sin
# argumentos se busca en output/results/entorno_parte1.RData y los resultados
# se escriben en output/.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. INSTALACIÓN/CARGA DE PAQUETES
# -----------------------------------------------------------------------------
# Se comprueban e instalan los paquetes necesarios si no están presentes.
# Nota: BiocManager es necesario para paquetes de Bioconductor.

required_cran <- c("caret", "glmnet", "randomForest", "e1071",
                   "gbm", "pROC", "ggplot2", "dplyr", "doParallel")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Carga de librerías
library(caret)          # Framework principal de entrenamiento de modelos
# NOTA: genefu y otros paquetes de Bioconductor cargados en la parte I pueden
# enmascarar train(), trainControl() y confusionMatrix(). Por eso a lo largo del
# script se llaman con el prefijo caret:: y no dependen del orden de carga.
library(glmnet)         # Modelo LASSO / Ridge / ElasticNet
library(randomForest)   # Modelo Random Forest
library(e1071)          # SVM (interfaz para kernlab / libsvm)
library(gbm)            # Gradient Boosting Machine
library(pROC)           # Cálculo y visualización de curvas ROC / AUC
library(ggplot2)        # Visualización avanzada
library(dplyr)          # Manipulación de datos
library(doParallel)     # Paralelización del entrenamiento (reduce tiempo)

# Paralelización: se usan todos los núcleos menos 1 para no saturar el sistema
n_cores <- max(1, detectCores() - 1)
registerDoParallel(cores = n_cores)
cat("Núcleos utilizados para paralelización:", n_cores, "\n")

# Semilla de aleatoriedad para reproducibilidad de todos los análisis
set.seed(42)

# -----------------------------------------------------------------------------
# DIRECTORIO DE SALIDA
# -----------------------------------------------------------------------------
# Los archivos generados se guardan en output/ (figuras en plots/, tablas y
# objetos en results/). Se admite otra carpeta como segundo argumento.
args <- commandArgs(trailingOnly = TRUE)

OUTPUT_DIR  <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "output"
PLOTS_DIR   <- file.path(OUTPUT_DIR, "plots")
RESULTS_DIR <- file.path(OUTPUT_DIR, "results")
dir.create(PLOTS_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Funciones auxiliares para construir rutas de salida
plot_path   <- function(f) file.path(PLOTS_DIR, f)
result_path <- function(f) file.path(RESULTS_DIR, f)
cat("Todos los resultados se guardarán en:", normalizePath(OUTPUT_DIR), "\n\n")


# -----------------------------------------------------------------------------
# 1. CARGA DEL ENTORNO PREVIO (RData de la Unidad 5)
# -----------------------------------------------------------------------------
# Se carga el entorno guardado por 01_mammaprint_labels_and_lasso.R, que
# contiene los objetos X_ml, risk, train_df, test_df, model, cor_mat y
# high_cor_idx. Reutilizarlos garantiza que la partición train/test es la
# misma y que los modelos son comparables con el LASSO de la parte I.
rdata_path <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  file.path(RESULTS_DIR, "entorno_parte1.RData")
}

if (!file.exists(rdata_path)) {
  stop(sprintf(paste0("No se encuentra el entorno de la parte I: '%s'.\n",
                      "  Ejecuta antes 01_mammaprint_labels_and_lasso.R o pasa la ruta como argumento."),
               rdata_path), call. = FALSE)
}

cat("Entorno de la parte I:", normalizePath(rdata_path), "\n")
load(rdata_path)

# Inspección de los objetos cargados en el entorno
cat("\n--- Objetos cargados desde el RData ---\n")
print(ls())

# ============================================================
# ASIGNACIÓN DE OBJETOS DEL ENTORNO PREVIO (Unidad 5)
# ============================================================
# Tras inspeccionar el entorno cargado, los objetos relevantes son:
#
#   X_ml     : matrix 75 x 1000
#              Matriz de expresión filtrada a los 1000 genes de mayor
#              varianza, ya normalizada y lista para ML. 
#
#   risk     : character 75
#              Etiquetas MammaPrint ("High Risk" / "Low Risk") para
#              cada una de las 75 muestras luminales.
#
#   train_df : data.frame 53 x 1001
#   test_df  : data.frame 22 x 1001
#              Partición previa (train/test) que incluye las 1000
#              columnas de genes + columna de etiqueta de clase.
#              Se reutilizan para mantener coherencia con la actividad anterior.
# ============================================================

# Matriz de predictores (muestras x genes)
expr_data <- X_ml
cat("\nDimensiones de X_ml (expr_data):", dim(expr_data), "\n")

# Vector de etiquetas, convertir a factor con nombres válidos para caret
# make.names() convierte "High Risk" -> "High.Risk", "Low Risk" -> "Low.Risk"
# Se convierten las clases a factor porque caret requiere variable respuesta
# categórica para clasificación supervisada binaria.
mp_labels <- factor(risk)
levels(mp_labels) <- make.names(levels(mp_labels))
cat("Etiquetas (risk) convertidas a factor:\n")
print(table(mp_labels))

# Verificación de alineación entre matriz y etiquetas
stopifnot(nrow(expr_data) == length(mp_labels))
cat("✔ Alineación matriz-etiquetas: CORRECTA\n")

# El objeto 'model' entrenado en la actividad anterior (glmnet) también está disponible
# y se incluirá en la comparación final de modelos como referencia
cat("Modelo previo Unidad 5 disponible como 'model':", class(model)[1], "\n")


# -----------------------------------------------------------------------------
# 2. PREPARACIÓN DE DATOS PARA CARET
# -----------------------------------------------------------------------------
# X_ml ya viene preprocesada:
#   - Filtrada a los 1000 genes de mayor varianza
#   - Alta correlación ya tratada (cor_mat y high_cor_idx están en el entorno)
#   - Normalizada con logCPM (edgeR)
#   - Muestras en filas, genes en columnas
#
# Solo se realizan los pasos mínimos necesarios para caret:
#   - Asegurar nombres de columna válidos en R
#   - Confirmar que no quedan NAs

cat("\n--- Preparación de datos para caret ---\n")
cat("Datos procedentes de X_ml (ya filtrados en Unidad 5)\n")
cat("Dimensiones originales:", dim(expr_data), "\n")

# Verificar ausencia de NAs
n_na <- sum(is.na(expr_data))
cat("Valores NA en la matriz:", n_na, "\n")
if (n_na > 0) stop("ERROR: hay NAs en la matriz de expresión. Revisar preprocesado.")

# Los diagnósticos de varianza casi nula y de correlación se realizaron en la
# parte I sobre el conjunto de entrenamiento; se documentan aquí para
# reproducibilidad. X_ml conserva los 1000 genes de mayor varianza.
cat("Diagnósticos realizados en la parte I:\n")
cat("  - Varianza casi nula evaluada con nearZeroVar (nzv_metrics disponible)\n")
cat("  - Genes con correlación > 0.90 identificados:", length(high_cor_idx), "genes\n")
cat("  - Genes finales en X_ml:", ncol(expr_data), "\n")

# -----------------------------------------------------------------------------
# 3. PARTICIONADO DE DATOS: ENTRENAMIENTO Y VALIDACIÓN
# -----------------------------------------------------------------------------
# La partición ya fue realizada en la actividad anterior y está guardada en el entorno
# como train_df (53 muestras) y test_df (22 muestras), ambos con 1000 genes
# más la columna de etiqueta de clase.
#
# Se reutiliza esta partición para mantener coherencia con el análisis previo
# y garantizar que el set de validación es completamente independiente.
# La proporción resultante es ~70/30, adecuada para el tamaño del dataset.
#
# NOTA: Reutilizar la misma partición evita el data leakage y asegura que
# los resultados son comparables con el modelo LASSO de la actividad anterior.

cat("\n--- Reutilizando partición train/test de la actividad anterior ---\n")

# La columna de etiquetas se llama "y", primera columna de train_df/test_df
# NOTA: se reutiliza la partición de la actividad anterior en lugar de generar una nueva
# para mantener coherencia analítica y evitar comparaciones no equivalentes.
label_col <- "y"  # La columna de etiquetas se llama "y" (primera columna de train_df/test_df)

# Separar predictores y etiquetas en entrenamiento
X_train <- as.data.frame(train_df[, setdiff(colnames(train_df), label_col)])
y_train <- factor(make.names(train_df[[label_col]]))

# Separar predictores y etiquetas en validación
X_val <- as.data.frame(test_df[, setdiff(colnames(test_df), label_col)])
y_val <- factor(make.names(test_df[[label_col]]), levels = levels(y_train))

# Verificacion de alineacion (detiene el script si hay desajuste)
stopifnot(nrow(X_train) == length(y_train))
stopifnot(nrow(X_val)   == length(y_val))
cat("Alineacion X_train/y_train:", nrow(X_train), "muestras\n")
cat("Alineacion X_val/y_val   :", nrow(X_val),   "muestras\n")

cat("Muestras en entrenamiento:", nrow(X_train), "\n")
print(table(y_train))
cat("Muestras en validación:", nrow(X_val), "\n")
print(table(y_val))


# -----------------------------------------------------------------------------
# 4. DEFINICIÓN DEL ESQUEMA DE VALIDACIÓN CRUZADA
# -----------------------------------------------------------------------------
# Se utiliza validación cruzada repetida (k=5, repeticiones=3) para estimar
# el rendimiento de los modelos durante el entrenamiento y optimizar los
# hiperparámetros de forma robusta, minimizando la varianza en la estimación.
#
# Además se activa classProbs = TRUE para poder calcular curvas ROC/AUC
# y summaryFunction = twoClassSummary para usar ROC como métrica principal.

# NOTA SOBRE ELECCIÓN DE FOLDS:
# Con n=53 muestras de entrenamiento, usar k=5 folds deja ~10 muestras por fold,
# lo que puede generar folds con una sola clase representada. Se reduce a k=3 folds
# (~17 muestras por fold) para mayor estabilidad, manteniendo 5 repeticiones
# para compensar la varianza de la estimación.
ctrl <- caret::trainControl(
  method          = "repeatedcv",   # Validación cruzada repetida
  number          = 3,              # 3 folds (más estable con n=53)
  repeats         = 5,              # 5 repeticiones para reducir varianza
  classProbs      = TRUE,           # Necesario para calcular AUC/ROC
  summaryFunction = twoClassSummary,# Métricas: ROC, Sens, Spec
  savePredictions = "final",        # Guardar predicciones del fold final
  verboseIter     = FALSE,          # Sin output por iteración
  allowParallel   = TRUE            # Aprovechar paralelización
)

cat("\n--- Control de entrenamiento configurado ---\n")
cat("Método: repeated CV | Folds:", 3, "| Repeticiones:", 5, "\n")
cat("Justificación: k=3 elegido por tamaño reducido del dataset (n=53 train)\n")
cat("Métrica de optimización: ROC (AUC)\n")


# =============================================================================
# 5. ENTRENAMIENTO DE MODELOS
# =============================================================================
# Se entrenan 4 modelos con distintos enfoques:
#   (A) LASSO     - Modelo lineal con regularización L1 (selección de genes)
#   (B) Random Forest - Ensemble de árboles, robusto y no lineal
#   (C) SVM lineal    - Margen máximo, eficaz en espacios de alta dimensión
#   (D) GBM           - Boosting de árboles, muy competitivo en clasificación


# -----------------------------------------------------------------------------
# 5A. LASSO (glmnet con alpha = 1)
# -----------------------------------------------------------------------------
# LASSO (Least Absolute Shrinkage and Selection Operator) aplica penalización
# L1 sobre los coeficientes. Muchos coeficientes se anulan-selección
# automática de genes relevantes. Ideal para datos de alta dimensionalidad
# como expresión génica.
#
# Hiperparámetro a optimizar:
#   - lambda: fuerza de regularización (valores más altos → más genes a 0)
#   - alpha = 1 fijado (LASSO puro; alpha = 0 sería Ridge)

cat("\n=== [1/4] Entrenando LASSO (glmnet) ===\n")

# Rejilla de hiperparámetros: se prueban 20 valores de lambda entre
# 10^-4 y 10^1 en escala logarítmica
lambda_grid <- expand.grid(
  alpha  = 1,                                # LASSO puro
  lambda = 10^seq(-4, 1, length.out = 20)   # 20 valores de lambda
)

set.seed(42)
model_lasso <- caret::train(
  x          = X_train,
  y          = y_train,
  method     = "glmnet",          # Regresión regularizada (LASSO/Ridge/EN)
  trControl  = ctrl,
  tuneGrid   = lambda_grid,
  metric     = "ROC",             # Optimizar por AUC
  preProcess = c("center", "scale")  # Estandarización: media 0, sd 1
  # centrar y escalar es obligatorio para LASSO para que la penalización
  # sea equitativa entre genes con distintas escalas de expresión
)

cat("✔ LASSO entrenado\n")
cat("  Mejor lambda:", model_lasso$bestTune$lambda, "\n")
cat("  AUC (CV):", max(model_lasso$results$ROC), "\n")

# Número de genes con coeficiente != 0 en el modelo óptimo
lasso_coefs <- coef(model_lasso$finalModel, s = model_lasso$bestTune$lambda)
n_genes_lasso <- sum(lasso_coefs != 0) - 1   # -1 para excluir intercepto
cat("  Genes seleccionados por LASSO:", n_genes_lasso, "\n")


# -----------------------------------------------------------------------------
# 5B. RANDOM FOREST
# -----------------------------------------------------------------------------
# Random Forest construye múltiples árboles de decisión sobre submuestras
# aleatorias de los datos (bagging) y en cada nodo usa un subconjunto
# aleatorio de predictores (feature randomness). La predicción final se
# obtiene por votación mayoritaria.
#
# Ventajas: robusto frente a overfitting, estima importancia de variables.
# Hiperparámetro a optimizar:
#   - mtry: número de variables candidatas en cada split del árbol
#           (por defecto: sqrt(p) para clasificación)

cat("\n=== [2/4] Entrenando Random Forest ===\n")

# Rejilla de mtry: se prueban 4 valores alrededor de sqrt(ncol(X_train))
p <- ncol(X_train)
mtry_grid <- expand.grid(
  mtry = c(floor(sqrt(p)/2), floor(sqrt(p)), floor(sqrt(p)*2), floor(sqrt(p)*3))
)

set.seed(42)
model_rf <- caret::train(
  x         = X_train,
  y         = y_train,
  method    = "rf",               # Random Forest (paquete randomForest)
  trControl = ctrl,
  tuneGrid  = mtry_grid,
  metric    = "ROC",
  ntree     = 500,                # Número de árboles (parámetro fijo, no optimizado)
  importance = TRUE               # Calcular importancia de variables
 
)

cat("✔ Random Forest entrenado\n")
cat("  Mejor mtry:", model_rf$bestTune$mtry, "\n")
cat("  AUC (CV):", max(model_rf$results$ROC), "\n")


# -----------------------------------------------------------------------------
# 5C. SVM LINEAL (Support Vector Machine)
# -----------------------------------------------------------------------------
# SVM busca el hiperplano de máximo margen que separa las dos clases.
# La versión lineal es adecuada cuando los datos son linealmente separables
# o casi, como suele ocurrir en espacios de alta dimensión (muchos genes).
#
# Hiperparámetro a optimizar:
#   - C (cost): controla el trade-off entre margen amplio y errores de
#     clasificación. C pequeño - margen amplio (más regularización);
#     C grande → menos errores en train (riesgo de overfitting)

cat("\n=== [3/4] Entrenando SVM Lineal ===\n")

svm_grid <- expand.grid(
  C = c(0.001, 0.01, 0.1, 1, 10)  # 5 valores de cost en escala log
)

set.seed(42)
model_svm <- caret::train(
  x          = X_train,
  y          = y_train,
  method     = "svmLinear",        # SVM con kernel lineal (paquete kernlab)
  trControl  = ctrl,
  tuneGrid   = svm_grid,
  metric     = "ROC",
  preProcess = c("center", "scale") # Estandarización: necesaria para SVM
  # NOTA: SVM es sensible a la escala de los predictores, por lo que
  # la estandarización es imprescindible para un rendimiento óptimo
)

cat("✔ SVM Lineal entrenado\n")
cat("  Mejor C:", model_svm$bestTune$C, "\n")
cat("  AUC (CV):", max(model_svm$results$ROC), "\n")


# -----------------------------------------------------------------------------
# 5D. GBM (Gradient Boosting Machine)
# -----------------------------------------------------------------------------
# GBM construye un ensemble de árboles de forma secuencial: cada árbol
# nuevo intenta corregir los errores del ensemble anterior. Es uno de los
# algoritmos más potentes para tablas de datos estructurados.
#
# Hiperparámetros a optimizar:
#   - n.trees          : número de árboles (iteraciones de boosting)
#   - interaction.depth: profundidad máxima de cada árbol
#   - shrinkage        : tasa de aprendizaje (learning rate)
#   - n.minobsinnode   : mínimo de obs. en nodo terminal (regularización)
#
# NOTA SOBRE REJILLA REDUCIDA:
# Con n=53 muestras en train y k=3 folds, cada fold tiene ~17 muestras.
# Una rejilla grande (18 combinaciones) genera muchos folds inestables donde
# alguna clase puede quedar sin representación, produciendo ROC=NA. Se reduce la rejilla:
#   - n.minobsinnode = 5 (en lugar de 10) para tolerar folds pequeños
#   - shrinkage fijo en 0.1 (learning rate más informativo con pocos datos)
#   - interaction.depth máximo de 2 (árboles poco profundos, menos overfitting)
# Esto es una decisión tomada ante la limitación del tamaño muestral.

cat("\n=== [4/4] Entrenando GBM (Gradient Boosting) ===\n")

gbm_grid <- expand.grid(
  n.trees           = c(50, 100, 200),   # Árboles: valores moderados
  interaction.depth = c(1, 2),           # Profundidad baja (datos pequeños)
  shrinkage         = 0.1,              # Learning rate único y estable
  n.minobsinnode    = 5                 # Reducido para tolerar folds pequeños
)
cat("Combinaciones en rejilla GBM:", nrow(gbm_grid), "\n")

# trainControl específico para GBM con allowParallel desactivado
# (gbm puede tener conflictos con algunos backends de paralelización)
ctrl_gbm <- caret::trainControl(
  method          = "repeatedcv",
  number          = 3,
  repeats         = 5,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  verboseIter     = FALSE,
  allowParallel   = FALSE            # GBM más estable sin paralelización
)

set.seed(42)
model_gbm <- caret::train(
  x         = X_train,
  y         = y_train,
  method    = "gbm",
  trControl = ctrl_gbm,
  tuneGrid  = gbm_grid,
  metric    = "ROC",
  verbose   = FALSE
)

cat("✔ GBM entrenado\n")
cat("  Mejor configuración:\n")
print(model_gbm$bestTune)
cat("  AUC (CV):", max(model_gbm$results$ROC, na.rm = TRUE), "\n")


# =============================================================================
# 6. COMPARACIÓN DE MODELOS
# =============================================================================
# Se comparan los 4 modelos mediante remuestreo sobre las métricas de
# validación cruzada (ROC, Sensitivity, Specificity).

cat("\n=== Comparación de modelos (métricas CV) ===\n")

model_list <- caret::resamples(list(
  LASSO = model_lasso,
  RF    = model_rf,
  SVM   = model_svm,
  GBM   = model_gbm
))

# Resumen estadístico de las métricas por modelo
summary_resamps <- summary(model_list)
print(summary_resamps)

# Plot de comparación de modelos (boxplot de ROC por modelo)
# NOTA: bwplot de caret no acepta xlab/ylab directamente porque internamente
# ya genera esos labels. Se usa scales para ajustar el eje y se añade
# el título mediante el argumento main propio de lattice.
cat("\n--- Generando plot de comparación de modelos ---\n")
png(plot_path("comparacion_modelos_boxplot.png"), width = 900, height = 600, res = 120)
bwplot(model_list,
       metric = "ROC",
       main   = "Comparación de AUC (Validación Cruzada Repetida 5x3)")
dev.off()
cat("✔ Guardado: comparacion_modelos_boxplot.png\n")

# Plot de medias con IC 95% por modelo
png(plot_path("comparacion_modelos_diff.png"), width = 900, height = 600, res = 120)
dotplot(model_list,
        metric = "ROC",
        main   = "AUC medio con IC 95% por modelo")
dev.off()
cat("✔ Guardado: comparacion_modelos_diff.png\n")



# Selección del mejor modelo basado en AUC medio en CV
auc_means <- summary_resamps$statistics$ROC[, "Mean"]
best_model_name <- names(which.max(auc_means))
cat("\n>>> Modelo seleccionado:", best_model_name,
    "| AUC CV:", round(max(auc_means), 4), "\n")

# Asignar el mejor modelo a una variable genérica para la evaluación final
best_model <- get(paste0("model_", tolower(best_model_name)))


# =============================================================================
# 7. EVALUACIÓN FINAL EN EL SET DE VALIDACIÓN INDEPENDIENTE
# =============================================================================
# El set de validación no ha participado en ningún paso del entrenamiento
# ni en la optimización de hiperparámetros. Permite estimar el rendimiento
# real (generalización) del modelo seleccionado.

cat("\n=== Evaluación del modelo seleccionado en set de VALIDACIÓN ===\n")

# 7.1 Predicciones de clase y probabilidades
pred_class <- predict(best_model, newdata = X_val)
pred_probs <- predict(best_model, newdata = X_val, type = "prob")

# Determinar cuál columna de probabilidades corresponde a la clase positiva
# (asumimos que el primer nivel del factor es la clase positiva)
positive_class <- levels(y_val)[1]
cat("Clase positiva para AUC:", positive_class, "\n")

# 7.2 Matriz de confusión con métricas derivadas
cat("\n--- Matriz de Confusión ---\n")
conf_mat <- caret::confusionMatrix(pred_class, y_val, positive = positive_class)
print(conf_mat)

# 7.3 Extracción de métricas 
accuracy    <- conf_mat$overall["Accuracy"]
kappa       <- conf_mat$overall["Kappa"]
sensitivity <- conf_mat$byClass["Sensitivity"]
specificity <- conf_mat$byClass["Specificity"]
ppv         <- conf_mat$byClass["Pos Pred Value"]
npv         <- conf_mat$byClass["Neg Pred Value"]
f1          <- conf_mat$byClass["F1"]

cat("\n--- Métricas de rendimiento ---\n")
cat(sprintf("  Accuracy    : %.4f\n", accuracy))
cat(sprintf("  Kappa       : %.4f\n", kappa))
cat(sprintf("  Sensitivity : %.4f\n", sensitivity))
cat(sprintf("  Specificity : %.4f\n", specificity))
cat(sprintf("  PPV         : %.4f\n", ppv))
cat(sprintf("  NPV         : %.4f\n", npv))
cat(sprintf("  F1-Score    : %.4f\n", f1))

# 7.4 Curva ROC y AUC en validación
roc_obj <- roc(response  = y_val,
               predictor = pred_probs[[positive_class]],
               levels    = rev(levels(y_val)))  # clase positiva vs negativa
auc_val <- auc(roc_obj)
cat(sprintf("\n  AUC (validación): %.4f\n", auc_val))
cat(sprintf("  IC 95%% AUC     : %s\n",
            paste(round(ci.auc(roc_obj), 4), collapse = " - ")))


# =============================================================================
# 8. VISUALIZACIONES FINALES
# =============================================================================

# 8.1 Curva ROC en set de validación
cat("\n--- Generando curva ROC ---\n")
png(plot_path("curva_ROC_validacion.png"), width = 700, height = 700, res = 120)
plot(roc_obj,
     col       = "#E63946",
     lwd       = 2.5,
     main      = paste0("Curva ROC - Modelo ", best_model_name,
                        " (Set de Validación)"),
     xlab      = "1 - Especificidad (Tasa Falsos Positivos)",
     ylab      = "Sensibilidad (Tasa Verdaderos Positivos)",
     print.auc = TRUE,
     auc.polygon = TRUE,
     auc.polygon.col = "#E6394620",
     grid      = TRUE)
abline(a = 0, b = 1, lty = 2, col = "grey50")   # Línea de azar
dev.off()
cat("✔ Guardado: curva_ROC_validacion.png\n")

# 8.2 Matriz de confusión como heatmap con ggplot2
cat("--- Generando heatmap de matriz de confusión ---\n")
conf_df <- as.data.frame(conf_mat$table)
colnames(conf_df) <- c("Predicho", "Real", "Freq")

p_conf <- ggplot(conf_df, aes(x = Real, y = Predicho, fill = Freq)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = Freq), size = 8, fontface = "bold", color = "white") +
  scale_fill_gradient(low = "#457B9D", high = "#E63946",
                      name = "N muestras") +
  labs(
    title    = paste0("Matriz de Confusión - ", best_model_name),
    subtitle = paste0("Set de validación (n = ", nrow(X_val), ")"),
    x        = "Clase Real",
    y        = "Clase Predicha",
    caption  = paste0("Accuracy = ", round(accuracy, 3),
                      " | AUC = ", round(auc_val, 3))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
    axis.text     = element_text(size = 12)
  )

ggsave(plot_path("matriz_confusion_heatmap.png"), p_conf,
       width = 7, height = 6, dpi = 150)
cat("✔ Guardado: matriz_confusion_heatmap.png\n")

# 8.3 Importancia de variables del mejor modelo
cat("--- Generando plot de importancia de variables ---\n")

# La importancia de variables se calcula de distinta forma según el modelo:
#   - RF  : MeanDecreaseGini / MeanDecreaseAccuracy
#   - GBM : Relative influence
#   - LASSO: valor absoluto de los coeficientes
#   - SVM : no tiene importancia nativa, se usa varImp genérico de caret

var_imp <- caret::varImp(best_model, scale = TRUE)  # Normalizado 0-100

# Representar top 20 variables más importantes
png(plot_path("importancia_variables.png"), width = 900, height = 700, res = 120)
plot(var_imp,
     top  = 20,
     main = paste0("Top 20 genes más importantes - Modelo ", best_model_name))
dev.off()
cat("✔ Guardado: importancia_variables.png\n")

# 8.4 (Extra) Curvas ROC superpuestas de los 4 modelos en validación
cat("--- Generando curvas ROC comparativas de todos los modelos ---\n")

# Función auxiliar para extraer probabilidades de cada modelo
get_roc <- function(model, X, y, pos_class) {
  probs <- predict(model, newdata = X, type = "prob")[[pos_class]]
  roc(response = y, predictor = probs, levels = rev(levels(y)))
}

roc_lasso <- get_roc(model_lasso, X_val, y_val, positive_class)
roc_rf    <- get_roc(model_rf,    X_val, y_val, positive_class)
roc_svm   <- get_roc(model_svm,   X_val, y_val, positive_class)
roc_gbm   <- get_roc(model_gbm,   X_val, y_val, positive_class)

png(plot_path("curvas_ROC_todos_modelos.png"), width = 800, height = 800, res = 120)
plot(roc_lasso, col = "#2196F3", lwd = 2,
     main = "Curvas ROC - Todos los Modelos (Set Validación)")
plot(roc_rf,    col = "#4CAF50", lwd = 2, add = TRUE)
plot(roc_svm,   col = "#FF9800", lwd = 2, add = TRUE)
plot(roc_gbm,   col = "#E63946", lwd = 2, add = TRUE)
abline(a = 0, b = 1, lty = 2, col = "grey60")
legend("bottomright",
       legend = c(
         paste0("LASSO (AUC=", round(auc(roc_lasso), 3), ")"),
         paste0("RF    (AUC=", round(auc(roc_rf),    3), ")"),
         paste0("SVM   (AUC=", round(auc(roc_svm),   3), ")"),
         paste0("GBM   (AUC=", round(auc(roc_gbm),   3), ")")
       ),
       col    = c("#2196F3", "#4CAF50", "#FF9800", "#E63946"),
       lwd    = 2,
       bty    = "n",
       cex    = 0.9)
dev.off()
cat("✔ Guardado: curvas_ROC_todos_modelos.png\n")


# =============================================================================
# 9. RESUMEN TABLA COMPARATIVA DE MÉTRICAS (VALIDACIÓN)
# =============================================================================
# Se construye una tabla con las métricas clave de todos los modelos en el
# set de validación.

cat("\n=== Tabla comparativa de métricas en validación ===\n")

evaluar_modelo <- function(model, X, y, pos_class) {
  pred_c  <- predict(model, newdata = X)
  pred_p  <- predict(model, newdata = X, type = "prob")[[pos_class]]
  cm      <- caret::confusionMatrix(pred_c, y, positive = pos_class)
  roc_obj <- roc(y, pred_p, levels = rev(levels(y)), quiet = TRUE)
  data.frame(
    Accuracy    = round(cm$overall["Accuracy"], 4),
    Kappa       = round(cm$overall["Kappa"],    4),
    Sensitivity = round(cm$byClass["Sensitivity"], 4),
    Specificity = round(cm$byClass["Specificity"], 4),
    F1          = round(cm$byClass["F1"],           4),
    AUC         = round(auc(roc_obj), 4),
    row.names   = NULL
  )
}

tabla_resultados <- rbind(
  data.frame(Modelo = "LASSO", evaluar_modelo(model_lasso, X_val, y_val, positive_class)),
  data.frame(Modelo = "RF",    evaluar_modelo(model_rf,    X_val, y_val, positive_class)),
  data.frame(Modelo = "SVM",   evaluar_modelo(model_svm,   X_val, y_val, positive_class)),
  data.frame(Modelo = "GBM",   evaluar_modelo(model_gbm,   X_val, y_val, positive_class))
)

print(tabla_resultados)

# Guardar tabla en CSV.
write.csv(tabla_resultados, result_path("tabla_comparativa_metricas.csv"),
          row.names = FALSE)
cat("✔ Guardado: tabla_comparativa_metricas.csv\n")


# =============================================================================
# 10. GUARDADO DEL ENTORNO Y OBJETOS CLAVE
# =============================================================================
# Se guardan todos los modelos y resultados para reproducibilidad

cat("\n--- Guardando entorno de la Unidad 6 ---\n")
save(
  model_lasso, model_rf, model_svm, model_gbm,
  best_model, best_model_name,
  tabla_resultados, conf_mat,
  roc_lasso, roc_rf, roc_svm, roc_gbm,
  X_train, X_val, y_train, y_val,
  file = result_path("UA6_MammaPrint_modelos.RData")
)
cat("✔ Guardado: UA6_MammaPrint_modelos.RData\n")


# =============================================================================
# 11. SESSION INFO (Reproducibilidad)
# =============================================================================
# Garantiza que cualquier persona pueda reproducir el análisis con las
# mismas versiones de R y paquetes.

cat("\n=== sessionInfo() ===\n")
sessionInfo()

# Detener clúster de paralelización al finalizar
stopImplicitCluster()

cat("\n=== ANÁLISIS COMPLETADO ===\n")
cat("Archivos generados en:", normalizePath(OUTPUT_DIR), "\n")
cat("  - comparacion_modelos_boxplot.png\n")
cat("  - comparacion_modelos_diff.png\n")
cat("  - curva_ROC_validacion.png\n")
cat("  - curvas_ROC_todos_modelos.png\n")
cat("  - matriz_confusion_heatmap.png\n")
cat("  - importancia_variables.png\n")
cat("  - tabla_comparativa_metricas.csv\n")
cat("  - UA6_MammaPrint_modelos.RData\n")
cat("\nCarpeta de resultados:", normalizePath(OUTPUT_DIR), "\n")

# =============================================================================
# FIN DEL SCRIPT
# =============================================================================
