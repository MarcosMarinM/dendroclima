# ==============================================================================
# 06_detrend_chron_build.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Comparative Chronology Builder.
#   Simultaneously applies three distinct detrending methods (Spline 30yr, 
#   Spline 100yr, and Modified Negative Exponential) to the same dataset and 
#   constructs a robust "Mean Median" chronology from the ensemble.
#
#   Methodology:
#   1. Detrending:
#      - Spline 30yr: High-frequency signal (annual variability).
#      - Spline 100yr: Medium-frequency signal (decadal trends).
#      - ModNegExp: Low-frequency/Biological trend removal.
#   2. Aggregation: Calculates the median of the three standardised indices for 
#      each year to create a consensus chronology less sensitive to method selection.
#
#   Inputs:
#   - RWL file.
#
#   Outputs:
#   - Individual Chronologies (TXT, PDF).
#   - Combined Median Chronology (TXT).
# ==============================================================================

# Cargar la librería dplR
library(dplR)

# --- PARÁMETROS A MODIFICAR ---
FILE_PATH_RWL <- "PLACEHOLDER/path/to/input_data.rwl"
OUTPUT_DIR    <- "PLACEHOLDER/path/to/output_dir"

# Cargar los datos desde el archivo
data <- read.rwl(FILE_PATH_RWL)

# Función para aplicar detrend a una serie y construir la cronología media
build_mean_chronology <- function(data, method, nyrs = NULL, prewhiten = FALSE) {
  detrended_data <- detrend(data, method = method, nyrs = nyrs)
  chronology <- chron(detrended_data, prewhiten = prewhiten)
  return(chronology)
}

# Construir cronologías residuales para cada método
chronology_spline30_residual <- build_mean_chronology(data, method = "Spline", nyrs = 30, prewhiten = TRUE)
chronology_spline100_residual <- build_mean_chronology(data, method = "Spline", nyrs = 100, prewhiten = TRUE)
chronology_modnegexp_residual <- build_mean_chronology(data, method = "ModNegExp", prewhiten = TRUE)

# Crear archivos PDF para guardar los gráficos
# Spline 30 años
pdf(file.path(OUTPUT_DIR, "cronologia_spline30.pdf"), width = 14, height = 8)
plot(chronology_spline30_residual, add.spline = TRUE, nyrs = 20, main = "Cronología residual (spline 30 años)", xlab = "Años", ylab = "RWI")
dev.off()

# Spline 100 años
pdf(file.path(OUTPUT_DIR, "cronologia_spline100.pdf"), width = 14, height = 8)
plot(chronology_spline100_residual, add.spline = TRUE, nyrs = 20, main = "Cronología residual (spline 100 años)", xlab = "Años", ylab = "RWI")
dev.off()

# ModNegExp
pdf(file.path(OUTPUT_DIR, "cronologia_modnegexp.pdf"), width = 14, height = 8)
plot(chronology_modnegexp_residual, add.spline = TRUE, nyrs = 20, main = "Cronología residual (ModNegExp)", xlab = "Años", ylab = "RWI")
dev.off()

# Función para guardar las cronologías como archivos .txt con la columna "year"
save_chronology <- function(chronology, file_path) {
  # Convertir los nombres de las filas en una columna llamada "year"
  chronology_df <- cbind(year = as.numeric(rownames(chronology)), as.data.frame(chronology))
  # Guardar el archivo
  write.table(chronology_df, file = file_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# Guardar los archivos
save_chronology(chronology_spline30_residual, file.path(OUTPUT_DIR, "cronologia_spline30.txt"))
save_chronology(chronology_spline100_residual, file.path(OUTPUT_DIR, "cronologia_spline100.txt"))
save_chronology(chronology_modnegexp_residual, file.path(OUTPUT_DIR, "cronologia_modnegexp.txt"))

# Leer las cronologías generadas
spline30 <- read.table(file.path(OUTPUT_DIR, "cronologia_spline30.txt"), header = TRUE, sep = "\t")
spline100 <- read.table(file.path(OUTPUT_DIR, "cronologia_spline100.txt"), header = TRUE, sep = "\t")
modnegexp <- read.table(file.path(OUTPUT_DIR, "cronologia_modnegexp.txt"), header = TRUE, sep = "\t")

# Combinar las columnas de interés (std, res) y mantener samp.depth
combine_cronologies <- function(spline30, spline100, modnegexp) {
  combined <- data.frame(
    year = spline30$year,  # Mantener los años
    std = apply(cbind(spline30$std, spline100$std, modnegexp$std), 1, function(x) median(x, na.rm = TRUE)),  # Mediana de std
    res = apply(cbind(spline30$res, spline100$res, modnegexp$res), 1, function(x) ifelse(all(is.na(x)), NA, median(x, na.rm = TRUE))),  # Mediana de res (ignorar NAs parciales)
    samp.depth = spline30$samp.depth  # Mantener igual
  )
  return(combined)
}

# Crear la cronología combinada
combined_cronology <- combine_cronologies(spline30, spline100, modnegexp)

# Guardar el archivo combinado
output_combined_path <- file.path(OUTPUT_DIR, "cronologia_mediana.txt")
write.table(combined_cronology, file = output_combined_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

