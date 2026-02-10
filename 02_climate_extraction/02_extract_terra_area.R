# ==============================================================================
# 03_extract_terra_area.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   TerraClimate Area Extractor (Polygon-based).
#   Extracts monthly climate data from TerraClimate NetCDF files for a defined 
#   geographic polygon (e.g., a sampling site, a watershed, or a protected area).
#
#   Methodology:
#   - Uses vector processing (sf/terra) to mask the climate grid with a user-defined 
#     polygon.
#   - Calculates the spatially weighted mean of the variable (e.g., precipitation, 
#     temperature) within that polygon for every month.
#   - Iterates through all NetCDF files in a directory to build a complete time series.
#
#   Inputs:
#   - TerraClimate NetCDF files.
#   - Polygon coordinates (Matrix).
#
#   Outputs:
#   - Time series text file (Monthly data).
# ==============================================================================

# --- PARÁMETROS A MODIFICAR ---
library(terra)
library(tidyverse)
library(sf)

# --- PARÁMETROS A MODIFICAR ---
CODIGO_VAR    <- "ppt" 
RUTA_ARCHIVOS <- "PLACEHOLDER/path/to/terraclimate_files/"
RUTA_SALIDA   <- "PLACEHOLDER/path/to/output_dir"

# Coordenadas del polígono (Matriz de puntos WGS84)
# Formato: c(lon, lat). Orden: Antihorario o Horario cerrando el polígono.
COORDS_POLIGONO <- matrix(c(
  -5.665032, 32.677370,
  -4.743459, 32.677370,
  -4.743459, 33.526580,
  -5.665032, 33.526580,
  -5.665032, 32.677370
), ncol = 2, byrow = TRUE)
# ------------------------------




# -----------------------------------------------------------------------------
# 3. Preparación del Entorno
# -----------------------------------------------------------------------------

# Detectar archivos .nc en la ruta especificada
patron_archivo <- paste0("\\.nc$")
archivos <- list.files(path = RUTA_ARCHIVOS, pattern = patron_archivo, full.names = TRUE)

# Generar nombre automático del archivo de salida
nombre_txt <- paste0("mediaarea_", CODIGO_VAR, "_terraclimate.txt")
archivo_final <- file.path(RUTA_SALIDA, nombre_txt)

print(paste0("--- INICIANDO PROCESO: ", toupper(CODIGO_VAR), " ---"))
print(paste("Archivos detectados:", length(archivos)))
print(paste("Salida prevista:", archivo_final))

# -----------------------------------------------------------------------------
# 4. Definición del Polígono (Coordenadas del área de estudio)
# -----------------------------------------------------------------------------
# Crear objeto vectorial sf/terra usando COORDS_POLIGONO
poligono <- st_polygon(list(COORDS_POLIGONO)) %>% 
  st_sfc(crs = 4326) %>% 
  vect()


# -----------------------------------------------------------------------------
# 5. Función de Procesamiento (Iterativa)
# -----------------------------------------------------------------------------

procesar_nc <- function(archivo) {
  
  # A. Extraer año del nombre del archivo (ej. "..._1958.nc")
  anio <- str_extract(basename(archivo), "\\d{4}")
  
  # B. Cargar Raster
  r <- rast(archivo)
  
  # Control de calidad básico: Verificar capas
  if(nlyr(r) > 12) { 
    warning(paste("El archivo", basename(archivo), "tiene más de 12 capas."))
  }
  
  # C. Extracción de datos
  # Usamos terra::extract con la media (fun = mean) ponderada implícita
  extraccion <- terra::extract(r, poligono, fun = mean, na.rm = TRUE)
  
  # D. Limpieza y formateo
  valores <- as.numeric(extraccion[ , -1]) # Eliminar columna ID
  fechas <- paste0(anio, sprintf("%02d", 1:12), "01") # Formato YYYYMM01
  
  # E. Gestión de memoria (Garbage Collection)
  rm(r)
  gc() 
  
  # F. Retorno
  return(data.frame(date = fechas, valor = valores))
}

# -----------------------------------------------------------------------------
# 6. Ejecución y Exportación
# -----------------------------------------------------------------------------

# Ejecutar bucle con manejo de errores
tryCatch({
  
  print("Procesando archivos... (Esto puede tardar unos segundos)")
  datos_completos <- map_dfr(archivos, procesar_nc)
  
  # Verificar si se generaron datos
  if(nrow(datos_completos) > 0) {
    
    # Ordenar cronológicamente
    datos_completos <- datos_completos %>% arrange(date)
    
    # Renombrar columna 'valor' al nombre real de la variable
    names(datos_completos)[2] <- CODIGO_VAR
    
    # Previsualización
    print("--- PRIMERAS FILAS ---")
    print(head(datos_completos))
    
    # Guardar a disco
    write.table(datos_completos, 
                file = archivo_final, 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
    
    print(paste("¡ÉXITO! Archivo guardado correctamente en:", archivo_final))
    
  } else {
    stop("La tabla de datos está vacía. Revisa las rutas.")
  }
  
}, error = function(e) {
  print(paste("ERROR CRÍTICO:", e$message))
})

# ==============================================================================
# FIN DEL SCRIPT
# ==============================================================================