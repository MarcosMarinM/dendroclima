# ==============================================================================
# 08_calc_spei.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Multi-Station SPEI Calculator (TerraClimate PET + Local Precip).
#   A batch processing script designed to calculate the Standardised Precipitation 
#   Evapotranspiration Index (SPEI) for a list of meteorological stations.
#
#   Workflow:
#   1. Hybrid Data Source: 
#      - Precipitation: Loaded from local station records (high accuracy).
#      - PET (Potential Evapotranspiration): Extracted from TerraClimate NetCDF 
#        files for the station coordinates (gap-free consistency).
#   2. Water Balance: Calculates P - PET for every month.
#   3. Multi-Scalar Analysis: Computes SPEI for all time scales from 1 to 48 months.
#
#   Inputs:
#   - Station Metadata (CSV with Lat/Lon).
#   - Station Precipitation Files.
#   - TerraClimate PET NetCDF folder.
#
#   Outputs:
#   - Individual Station SPEI files (scales 1-48).
# ==============================================================================

# 1. LIBRERÍAS
library(terra)
library(dplyr)
library(readr)
library(SPEI)

# --- PARÁMETROS A MODIFICAR ---
PATH_TERRACLIMATE   <- "PLACEHOLDER/path/to/terraclimate_pet/"
PATH_PRECIPITACION  <- "PLACEHOLDER/path/to/stations_precip/"
PATH_CSV_ESTACIONES <- "PLACEHOLDER/path/to/stations.csv"
PATH_SALIDA_SPEI    <- "PLACEHOLDER/path/to/output_spei/"
YEARS_ANALISIS      <- 1961:2019
# ------------------------------

if (!dir.exists(PATH_SALIDA_SPEI)) dir.create(PATH_SALIDA_SPEI, recursive = TRUE)

# -----------------------------------------------------------------------------
# 3. CARGA DE DATOS ESPACIALES (PET)
# -----------------------------------------------------------------------------
message(">>> Cargando metadatos y extrayendo PET...")

estaciones <- read_csv(PATH_CSV_ESTACIONES, show_col_types = FALSE)
vect_estaciones <- vect(estaciones, geom = c("lon", "lat"), crs = "EPSG:4326")

pet_acumulado <- vector("list", length(estaciones$est_id))
names(pet_acumulado) <- estaciones$est_id

years <- YEARS_ANALISIS 


# -- Bucle de extracción de PET (TerraClimate) --
for (yr in years) {
  archivo_nc <- file.path(PATH_TERRACLIMATE, paste0("TerraClimate_pet_", yr, ".nc"))
  
  if (file.exists(archivo_nc)) {
    r <- rast(archivo_nc)
    valores_meses <- terra::extract(r, vect_estaciones, ID = FALSE)
    
    for (i in 1:nrow(estaciones)) {
      id_est <- estaciones$est_id[i]
      pet_acumulado[[id_est]] <- c(pet_acumulado[[id_est]], as.numeric(valores_meses[i, ]))
    }
    cat("   - Procesado año:", yr, "\n")
  } else {
    warning(paste("   [!] Falta archivo:", archivo_nc))
    for (id_est in names(pet_acumulado)) {
      pet_acumulado[[id_est]] <- c(pet_acumulado[[id_est]], rep(NA, 12))
    }
  }
}

# -----------------------------------------------------------------------------
# 4. CÁLCULO DE SPEI MULTI-ESCALA (1-48)
# -----------------------------------------------------------------------------

message("\n>>> Calculando SPEI (escalas 1-48) por estación...")

escalas <- 1:48

for (i in 1:nrow(estaciones)) {
  
  id_actual <- estaciones$est_id[i]
  
  # Cargar Precipitación
  archivo_prec <- file.path(PATH_PRECIPITACION, paste0(id_actual, "_prec_estacion.txt"))
  
  if (file.exists(archivo_prec)) {
    df_prec <- read.table(archivo_prec, header = TRUE, sep = "\t")
    
    # Datos base
    P <- df_prec$precipitation
    PE <- pet_acumulado[[id_actual]]
    
    if (length(P) == length(PE)) {
      
      # Calcular Balance Hídrico (D = P - PET)
      balance <- P - PE
      
      # Crear matriz para guardar resultados
      # Filas = meses, Columnas = escalas (SPEI_1, SPEI_2...)
      matriz_spei <- matrix(NA, nrow = length(balance), ncol = length(escalas))
      colnames(matriz_spei) <- paste0("SPEI_", escalas)
      
      # --- BUCLE DE CÁLCULO POR ESCALA ---
      for (k in escalas) {
        tryCatch({
          # Calculamos SPEI para la escala k
          # verbose=FALSE para no llenar la consola
          spei_obj <- spei(balance, scale = k, verbose = FALSE)
          
          # Extraemos la serie temporal ajustada (fitted values)
          matriz_spei[, k] <- as.numeric(spei_obj$fitted)
          
        }, error = function(e) {
          message(paste("Error en estación", id_actual, "escala", k))
        })
      }
      
      # --- FORMATEO FINAL ---
      # Unir fecha con la matriz de resultados
      df_final <- data.frame(date = df_prec$date, matriz_spei)
      
      # Guardar archivo
      nombre_archivo <- file.path(PATH_SALIDA_SPEI, paste0(id_actual, "_SPEI_1_48.txt"))
      
      write.table(df_final, file = nombre_archivo, 
                  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
      
      cat("   [OK] Guardado:", id_actual, "\n")
      
    } else {
      warning(paste("Longitudes no coinciden para", id_actual))
    }
    
  } else {
    warning(paste("No se encuentra archivo de precipitación para:", id_actual))
  }
}

message("\n==============================================================================")
message(" CÁLCULO DE SPEI FINALIZADO ")
message(" Archivos guardados en: ", PATH_SALIDA_SPEI)
message("==============================================================================")
