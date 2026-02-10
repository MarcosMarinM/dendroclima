# ==============================================================================
# 04_detrend_adaptive.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Adaptive Detrending Algorithm.
#   Implements a logic-based standardisation approach where the detrending method 
#   is selected dynamically for each individual tree-ring series based on its 
#   biological age (length).
#
#   Methodology:
#   1. Classification by Age:
#      - Young (<50 years): Spline (50% length).
#      - Mature (50-100 years): Spline (67% length).
#      - Old (>100 years): Modified Negative Exponential (ModNegExp) to preserve 
#        long-term trends (geometric growth).
#   2. Spline Stiffness: Calculates 'nyrs' (spline stiffness) as a percentage 
#      of series length to avoid overfitting young trees or filtering out meaningful 
#      low-frequency signals in older trees.
#
#   Inputs: 
#   - Segmented RWL file (e.g., 'IBE_TRW_mid.rwl').
#
#   Outputs:
#   - Standardised Chronology (TXT, PDF).
#   - RWI file (Tucson).
#   - Signal Strength Statistics (SSS) plot.
# ==============================================================================

# Cargar la librería dplR
library(dplR)

# --- PARÁMETROS A MODIFICAR ---
FILE_PATH_RWL  <- 'PLACEHOLDER/path/to/input_data.rwl'
OUTPUT_DIR     <- 'PLACEHOLDER/path/to/output_dir'

# Nombres de salida (opcional modificar)
OUTPUT_TXT_NAME <- "output_file.txt"
OUTPUT_PDF_NAME <- "output_file.pdf"
OUTPUT_RWI_NAME <- "output_file.rwi"

# Cargar los datos desde el archivo
data <- read.rwl(FILE_PATH_RWL)

# Extraer los años reales del archivo .rwl
real_years <- as.numeric(rownames(data))

# Función para determinar el método de detrending en función de la edad de la serie
choose_detrend_method <- function(series_age) {
  if (series_age < 50) {
    return("Spline")  # Usar spline para series jóvenes
  } else if (series_age >= 50 & series_age < 100) {
    return("Spline")  # Usar spline para series de edad media
  } else {
    return("ModNegExp")  # Usar ModNegExp para series más viejas
  }
}

# Función para calcular nyrs adaptativo basado en la longitud de la serie
calculate_nyrs_adaptive <- function(series_length) {
  # Ajustar nyrs como porcentaje de la longitud de la serie:
  # - 67% para series > 100 años (equilibrio señal/ruido)
  # - 50% para series <= 100 años (evitar overfitting)
  # - Mínimo absoluto de 20 años (para series muy cortas)
  nyrs_custom <- ifelse(series_length > 100, 
                        round(0.67 * series_length), 
                        round(0.50 * series_length))
  
  nyrs_custom <- max(nyrs_custom, 20)  # Nunca menos de 20 años
  return(nyrs_custom)
}

# Función para aplicar detrend a una serie individual
detrend_individual_series <- function(series, method) {
  series_length <- length(series[!is.na(series)])  # Longitud de la serie sin NAs
  
  if (method == "Spline") {
    nyrs_custom <- calculate_nyrs_adaptive(series_length)
    detrended_series <- detrend.series(series, method = "Spline", nyrs = nyrs_custom)
  } else if (method == "ModNegExp") {
    detrended_series <- detrend.series(series, method = "ModNegExp")
  } else {
    stop("Método de detrending no reconocido")
  }
  return(detrended_series)
}

# Aplicar detrending a cada serie individualmente
detrended_series_list <- lapply(data, function(series) {
  series_age <- length(series[!is.na(series)])  # Longitud de la serie sin NAs
  method <- choose_detrend_method(series_age)
  detrended_series <- detrend_individual_series(series, method)
  return(detrended_series)
})

# Convertir la lista de series detrendizadas en un data frame
detrended_data <- as.data.frame(detrended_series_list)

# Asignar los años reales como nombres de fila
rownames(detrended_data) <- real_years

# Construir la cronología media
chronology <- chron(detrended_data, prewhiten = TRUE)

# Guardar la cronología como archivo .txt
save_chronology <- function(chronology, file_path) {
  # Convertir los nombres de las filas en una columna llamada "year"
  chronology_df <- cbind(year = as.numeric(rownames(chronology)), as.data.frame(chronology))
  # Guardar el archivo
  write.table(chronology_df, file = file_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

output_chronology_path <- file.path(OUTPUT_DIR, OUTPUT_TXT_NAME)
save_chronology(chronology, output_chronology_path)

# Crear un archivo PDF para guardar el gráfico de la cronología final
pdf(file.path(OUTPUT_DIR, OUTPUT_PDF_NAME), width = 14, height = 8)
plot(chronology, add.spline = TRUE, nyrs = 20, main = "Cronología final (detrending adaptativo)", xlab = "Años", ylab = "RWI")
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
output_rwi_path <- file.path(OUTPUT_DIR, OUTPUT_RWI_NAME)
save_rwi(detrended_data, output_rwi_path)

#--------------## ANÁLISIS DE SENSIBILIDAD DE LA SEÑAL (SSS) ##--------------#

# Asegúrate de que detrended_data tenga los años como nombres de fila
rownames(detrended_data) <- as.numeric(rownames(detrended_data))

# Aplicar el análisis SSS
sss_results <- sss(detrended_data)

# Mostrar los resultados del análisis SSS
print(sss_results)

# Graficar los resultados del análisis SSS con los años correctos en el eje x
plot(rownames(detrended_data), sss_results, type = "l", main = "Análisis de sensibilidad de la señal (SSS)", xlab = "Años", ylab = "SSS")
