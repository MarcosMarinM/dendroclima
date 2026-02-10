# ==============================================================================
# 09_calc_spei_area.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Area-Based SPEI Calculator.
#   Calculates the Standardised Precipitation Evapotranspiration Index (SPEI) 
#   representing a specific geographic region (Polygon), rather than a single point.
#
#   Methodology:
#   1. Input Synchronisation: Merges a pre-calculated regional Precipitation series 
#      with a regional PET series extracted from TerraClimate.
#   2. Water Balance: Computes P - PET.
#   3. SPEI Calculation: Generates SPEI for scales 1-48 months.
#
#   Inputs:
#   - Regional Precipitation File.
#   - TerraClimate PET NetCDF folder (for on-the-fly areal extraction).
#   - Polygon Coordinates.
#
#   Outputs:
#   - Regional SPEI file (scales 1-48).
# ==============================================================================

# -----------------------------------------------------------------------------
# 1. LIBRERÍAS
# -----------------------------------------------------------------------------
library(terra)
library(tidyverse)
library(SPEI)
library(sf)

# -----------------------------------------------------------------------------
# 2. CONFIGURACIÓN (¡MODIFICAR ESTO!)
# -----------------------------------------------------------------------------

# A) Archivo de Precipitación (Debe tener encabezado y al menos 2 columnas: date, value)
ARCHIVO_PRECIPITACION <- "PLACEHOLDER/path/to/precip_area.txt"

# B) Carpeta donde están los NetCDF de PET de TerraClimate (TerraClimate_pet_XXXX.nc)
DIR_PET_TERRACLIMATE <- "PLACEHOLDER/path/to/terraclimate_pet/"

# C) Ruta de salida para el archivo SPEI
RUTA_SALIDA <- "PLACEHOLDER/path/to/output_spei/"
NOMBRE_SALIDA <- "SPEI_Area_tc_1_48.txt"

if (!dir.exists(RUTA_SALIDA)) dir.create(RUTA_SALIDA, recursive = TRUE)

# D) Coordenadas del Polígono para extraer la PET (Mismo que usaste antes)
coords <- matrix(c(
  -4.92526842693917,  33.085650613862896,
  -4.884476923781775, 32.362264545687026,
  -5.646143799842935, 32.370080915391426,
  -5.498901634926672, 33.13259724802919,
  -4.92526842693917,  33.085650613862896
), ncol = 2, byrow = TRUE)

# -----------------------------------------------------------------------------
# 3. PROCESAMIENTO
# -----------------------------------------------------------------------------

# --- 3.1 Cargar Precipitación ---
message(">>> Cargando archivo de precipitación...")
df_prec <- read.table(ARCHIVO_PRECIPITACION, header = TRUE, sep = "\t")

# Asegurar nombres estándar (asumimos columna 1 es fecha, columna 2 es valor)
names(df_prec)[1] <- "date"
names(df_prec)[2] <- "P"
# Asegurar que date sea string para el join posterior (YYYYMM01)
df_prec$date <- as.character(df_prec$date) 

# --- 3.2 Extraer PET de TerraClimate (Media Areal) ---
message(">>> Iniciando extracción de PET areal (TerraClimate)...")

# Crear polígono
poligono <- st_polygon(list(coords)) %>% st_sfc(crs = 4326) %>% vect()

# Listar archivos PET
files_pet <- list.files(DIR_PET_TERRACLIMATE, pattern = "\\.nc$", full.names = TRUE)

# Función de extracción rápida (reutilizando lógica anterior)
extraer_pet_anio <- function(archivo) {
  anio <- str_extract(basename(archivo), "\\d{4}")
  r <- rast(archivo)
  # Extraer media ponderada del área
  val <- terra::extract(r, poligono, fun = mean, na.rm = TRUE)
  # Limpieza
  v_num <- as.numeric(val[, -1])
  fechas <- paste0(anio, sprintf("%02d", 1:12), "01")
  
  rm(r); gc() # Limpieza de memoria
  return(data.frame(date = fechas, PET = v_num))
}

# Ejecutar extracción en bloque
df_pet <- map_dfr(files_pet, extraer_pet_anio)

# --- 3.3 Fusión de Datos (Sync) ---
message(">>> Sincronizando precipitación y PET...")

# Inner Join: Solo se queda con los meses donde existen AMBOS datos.
# (Ej: Si Prec es 1901-2024 y PET es 1958-2020, se queda con 1958-2020)
df_clima <- inner_join(df_prec, df_pet, by = "date")

message(paste("   Periodo común encontrado:", nrow(df_clima), "meses."))
message(paste("   Desde:", min(df_clima$date), "Hasta:", max(df_clima$date)))

# -----------------------------------------------------------------------------
# 4. CÁLCULO DE SPEI (Escalas 1-48)
# -----------------------------------------------------------------------------

# Calcular Balance Hídrico (D = P - PET)
balance <- df_clima$P - df_clima$PET

# Matriz para resultados
escalas <- 1:48
matriz_spei <- matrix(NA, nrow = length(balance), ncol = length(escalas))
colnames(matriz_spei) <- paste0("SPEI_", escalas)

message(">>> Calculando índices SPEI...")

for (k in escalas) {
  tryCatch({
    # verbose = FALSE para limpiar consola
    # na.rm = TRUE permite calcular aunque falte algún dato aislado
    spei_calc <- spei(balance, scale = k, verbose = FALSE, na.rm = TRUE)
    
    # Extraemos la serie ajustada
    val_spei <- as.numeric(spei_calc$fitted)
    
    # IMPORTANTE: SPEI devuelve -Inf o Inf si hay ceros absolutos en distribuciones raras,
    # lo corregimos a NA o valores límite si fuera necesario, 
    # pero por defecto 'fitted' suele estar bien.
    matriz_spei[, k] <- val_spei
    
  }, error = function(e) {
    message(paste("   [!] Error calculando escala", k, ":", e$message))
  })
}

# -----------------------------------------------------------------------------
# 5. EXPORTACIÓN
# -----------------------------------------------------------------------------

# Unimos fechas con los valores de SPEI
df_final <- data.frame(date = df_clima$date, matriz_spei)

archivo_out_full <- file.path(RUTA_SALIDA, NOMBRE_SALIDA)

write.table(df_final, 
            file = archivo_out_full, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE, 
            na = "NA") # Los NA se escriben como texto "NA"

message("\n==============================================================================")
message(" ¡PROCESO COMPLETADO! ")
message(" Archivo guardado en: ", archivo_out_full)
message(" Dimensiones: ", nrow(df_final), " meses x ", ncol(df_final), " columnas")
message("==============================================================================")
