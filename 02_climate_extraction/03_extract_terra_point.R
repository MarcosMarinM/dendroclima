# ==============================================================================
# 04_extract_terra_point.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   TerraClimate Point Extractor.
#   Extracts monthly climate time series from TerraClimate NetCDF files for a 
#   single specific geographic coordinate (Pixel-based extraction).
#
#   Use Case:
#   - Ideal for pinpointing climate conditions for a specific tree or small stand 
#     where spatial averaging is not desired.
#   - Faster than areal extraction for single points.
#
#   Inputs:
#   - Directory of TerraClimate NetCDF files.
#   - Target Latitude/Longitude.
#
#   Outputs:
#   - Time series text file (Monthly data).
# ==============================================================================

# --- PARÁMETROS A MODIFICAR ---
library(terra)

# --- PARÁMETROS A MODIFICAR ---
RUTA_CARPETA   <- "PLACEHOLDER/path/to/terraclimate_files/"
ARCHIVO_SALIDA <- 'PLACEHOLDER/path/to/output_file.txt'

# Coordenadas del punto (WGS84)
MI_LON <- -5.07259
MI_LAT <- 33.0097
# ------------------------------


# --- 2. PROCESAMIENTO ITERATIVO ---

# Listar todos los archivos .nc
archivos_nc <- list.files(path = RUTA_CARPETA, pattern = "\\.nc$", full.names = TRUE)

# Ordenarlos por nombre para asegurar cronología (1958, 1959...)
archivos_nc <- sort(archivos_nc)

# Crear un contenedor para los resultados
# Cambiamos el nombre de la columna a 'vpd'
resultados_df <- data.frame(date = character(), vpd = numeric(), stringsAsFactors = FALSE)

cat("Iniciando extracción de VPD Procesando", length(archivos_nc), "archivos...\n")

# Bucle archivo por archivo
for (archivo in archivos_nc) {
  
  # 1. Cargar el raster
  r <- rast(archivo)
  
  # 2. Extraer el valor del punto (CORREGIDO AQUÍ)
  # Usamos terra::extract para evitar conflictos con otros paquetes
  valores <- terra::extract(r, cbind(MI_LON, MI_LAT))
  
  # 'valores' será un data.frame con 1 fila y 12 columnas
  vpd_vals <- as.numeric(t(valores))
  
  # 3. Generar las fechas
  year_str <- regmatches(basename(archivo), regexpr("\\d{4}", basename(archivo)))
  
  if (length(year_str) == 0) {
    warning(paste("No se pudo detectar el año en el archivo:", basename(archivo)))
    next
  }
  
  fechas <- sprintf("%s%02d01", year_str, 1:12)
  
  # 4. Unir y guardar
  temp_df <- data.frame(date = fechas, vpd = vpd_vals)
  resultados_df <- rbind(resultados_df, temp_df)
  
  rm(r, valores, vpd_vals)
}

# --- 3. GUARDAR RESULTADO FINAL ---

# Verificar si hay NAs
if (any(is.na(resultados_df$vpd))) {
  cat("¡Atención! Hay valores NA. Verifica que las coordenadas caigan en tierra firme.\n")
}

# Escribir en formato TXT separado por tabuladores
write.table(resultados_df, 
            file = ARCHIVO_SALIDA, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

cat("¡Proceso terminado! Datos guardados en:", ARCHIVO_SALIDA, "\n")
head(resultados_df) # Mostrar las primeras filas para verificar