# ==============================================================================
# 04_detrend_spline32.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   High-Frequency Standardisation (32-Year Spline).
#   Applies a fixed 32-year cubic smoothing spline to all series.
#
#   Purpose:
#   - To act as a high-pass filter, removing all low-frequency trends (biological 
#     growth, decadal climate trends) and retaining only inter-annual to 
#     short-term variability.
#   - Ideal for analyzing high-frequency climate responses (e.g., extreme events, 
#     annual correlation).
#   - Corresponds to a standard dendrochronological "high-pass" configuration.
#
#   Inputs: 
#   - RWL file.
#
#   Outputs:
#   - 32-year Spline Standardised Chronology (TXT, PDF).
#   - RWI file.
# ==============================================================================

# Cargar la librería dplR
library(dplR)

# --- PARÁMETROS A MODIFICAR ---
FILE_PATH_RWL <- 'PLACEHOLDER/path/to/input_data.rwl'
OUTPUT_DIR    <- 'PLACEHOLDER/path/to/output_dir'

# Cargar los datos desde el archivo
data <- read.rwl(FILE_PATH_RWL)

# Extraer los años reales del archivo .rwl
real_years <- as.numeric(rownames(data))

# Función para aplicar detrending adaptativo según la longitud de la serie
detrend_adaptativo <- function(series) {
  series_length <- length(series[!is.na(series)])  # Longitud de la serie sin NAs
  
  if (series_length >= 100) {
    detrended_series <- detrend.series(series, method = "Spline", nyrs = 32)
  } else if (series_length >= 50) {
    detrended_series <- detrend.series(series, method = "Spline", nyrs = 20)
  } else if (series_length >= 30) {
    # Intentar con ModNegExp, si no converge, usar spline de 10 años
    detrended_series <- tryCatch(
      detrend.series(series, method = "ModNegExp"),
      error = function(e) detrend.series(series, method = "Spline", nyrs = 10)
    )
  } else {
    # Series muy cortas: intentar ModNegExp, y si falla, usar un spline corto
    detrended_series <- tryCatch(
      detrend.series(series, method = "ModNegExp"),
      error = function(e) detrend.series(series, method = "Spline", nyrs = 5)
    )
  }
  
  return(detrended_series)
}

# Aplicar detrending a cada serie individualmente
detrended_series_list <- lapply(data, detrend_adaptativo)

# Convertir a data frame y mantener los años reales como nombres de fila
detrended_data <- as.data.frame(detrended_series_list)
rownames(detrended_data) <- real_years

# Construir la cronología media
chronology_modnegexp <- chron(detrended_data, prewhiten = TRUE)

# Guardar la cronología como archivo .txt
save_chronology <- function(chronology, file_path) {
  # Convertir los nombres de las filas en una columna llamada "year"
  chronology_df <- cbind(year = as.numeric(rownames(chronology)), as.data.frame(chronology))
  # Guardar el archivo
  write.table(chronology_df, file = file_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

output_chronology_path <- file.path(OUTPUT_DIR, output_chronology_txt)
save_chronology(chronology_modnegexp, output_chronology_path)

# Crear un archivo PDF para guardar el gráfico de la cronología final
pdf(file.path(OUTPUT_DIR, output_chronology_pdf), width = 14, height = 8)
plot(chronology_modnegexp, add.spline = TRUE, nyrs = 20, main = "Cronología final", xlab = "Años", ylab = "RWI")
dev.off()

#--------------## FUNCIÓN PARA GUARDAR EL .RWI ##--------------#
save_rwi <- function(detrended_data, file_path) {
  # Asegurar que los nombres de las columnas sean correctos (sin prefijos)
  colnames(detrended_data) <- gsub("\\.dat$", "", colnames(detrended_data))
  
  # Convertir a objeto "rwl" para compatibilidad con dplR
  detrended_rwl <- as.rwl(detrended_data)
  
  # Guardar como archivo .rwl (que técnicamente es un .rwi)
  write.rwl(detrended_rwl, fname = file_path, format = "tucson")
}

#--------------## USO DE LA FUNCIÓN ##--------------#
output_rwi_path <- file.path(OUTPUT_DIR, output_rwi)
save_rwi(detrended_data, output_rwi_path)

#--------------## ANÁLISIS DE SENSIBILIDAD DE LA SEÑAL (SSS) ##--------------#

# Asegúrate de que detrended_data tenga los años como nombres de fila
rownames(detrended_data) <- as.numeric(rownames(detrended_data))

# Aplicar el análisis SSS
sss_results <- sss(detrended_data)

# Mostrar los resultados del análisis SSS
print(sss_results)

# Graficar los resultados del análisis SSS con los años correctos en el eje x
pdf(file.path(OUTPUT_DIR, output_sss_pdf), width = 14, height = 8)
plot(rownames(detrended_data), sss_results, type = "l", main = "Análisis de sensibilidad de la señal", xlab = "Años", ylab = "SSS")
dev.off()
