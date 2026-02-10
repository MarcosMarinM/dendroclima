# ==============================================================================
# 11_calc_spi.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Standardised Precipitation Index (SPI) Calculator.
#   A batch processing script designed to calculate the SPI (Standardised Precipitation 
#   Index) for a list of meteorological stations using only precipitation data.
#   Similar to SPEI but relying solely on precipitation inputs.
#
#   Methodology:
#   - Fits a probability distribution (typically Gamma) to the precipitation series.
#   - Transforms probabilities to a standard normal distribution (Z-score).
#   - Calculates indices for multiple time scales (1-48 months).
#
#   Inputs:
#   - Directory of Station Precipitation Files.
#
#   Outputs:
#   - Individual Station SPI files (scales 1-48).
# ==============================================================================

# 1. LIBRERÍAS
library(dplyr)
library(readr)
library(SPEI)

# --- PARÁMETROS A MODIFICAR ---
PATH_PRECIPITACION <- "PLACEHOLDER/path/to/stations_precip/"
PATH_SALIDA_SPI    <- "PLACEHOLDER/path/to/output_spi/"
# ------------------------------ 

# Crear directorio si no existe
if (!dir.exists(PATH_SALIDA_SPI)) dir.create(PATH_SALIDA_SPI, recursive = TRUE)

# 3. PROCESAMIENTO
# Obtenemos la lista de archivos de estaciones
archivos_estaciones <- list.files(PATH_PRECIPITACION, pattern = "_prec_estacion.txt", full.names = TRUE)

cat(paste(">>> Encontrados", length(archivos_estaciones), "archivos de estación.\n"))
cat(">>> Iniciando cálculo de SPI (escalas 1-48)...\n\n")

# Definir escalas
escalas <- 1:48

for (archivo in archivos_estaciones) {
  
  # Extraer ID de la estación
  nombre_archivo <- basename(archivo)
  id_estacion <- sub("_prec_estacion.txt", "", nombre_archivo)
  
  # Leer datos
  df_prec <- read.table(archivo, header = TRUE, sep = "\t")
  P <- df_prec$precipitation
  
  # Matriz para guardar resultados
  matriz_spi <- matrix(NA, nrow = length(P), ncol = length(escalas))
  colnames(matriz_spi) <- paste0("SPI_", escalas)
  
  # --- BUCLE DE CÁLCULO (SPI) ---
  for (k in escalas) {
    tryCatch({
      # Cálculo SPI (Gamma default)
      spi_obj <- spi(P, scale = k, verbose = FALSE, distribution = "Gamma")
      matriz_spi[, k] <- as.numeric(spi_obj$fitted)
      
    }, error = function(e) {
      message(paste("   [Error] Estación", id_estacion, "- Escala", k, ":", e$message))
    })
  }
  
  # --- GUARDADO ---
  df_final <- data.frame(date = df_prec$date, matriz_spi)
  
  nombre_salida <- file.path(PATH_SALIDA_SPI, paste0(id_estacion, "_SPI_1_48.txt"))
  
  write.table(df_final, file = nombre_salida, 
              sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  
  cat("   [OK] Procesada:", id_estacion, "\n")
}

message("\n==============================================================================")
message(" CÁLCULO DE SPI FINALIZADO ")
message(" Archivos guardados en: ", PATH_SALIDA_SPI)
message("==============================================================================")
