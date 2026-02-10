# ==============================================================================
# 01_extract_phyda.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Paleo Corrected Hydroclimate Data Atlas (PHYDA) Extractor.
#   Designed to extract hydroclimate reconstruction data (PDSI, SPEI, Temperature) 
#   from the massive PHYDA NetCDF dataset for a specific region of interest.
#
#   Capabilities:
#   - Spatial Averaging: Calculates the mean value of climate variables within a 
#     defined lat/lon box.
#   - Nearest Neighbor Fallback: If the defined box is smaller than a grid cell, 
#     it intelligently selects the nearest pixel.
#   - Coordinate Handling: Automatically detects and adjusts for 0-360 vs -180/180 
#     longitude formats.
#   - Global Indices: Also extracts global climate indices (Nino, AMO, PDO) stored 
#     within the PHYDA dataset.
#
#   Inputs:
#   - PHYDA NetCDF file.
#   - Target Coordinates (Lat/Lon Box).
#
#   Outputs:
#   - Tab-separated text file with annual series for the selected region.
# ==============================================================================

# --- PARÁMETROS A MODIFICAR ---
library(ncdf4)

# --- PARÁMETROS A MODIFICAR ---
NC_PATH     <- 'PLACEHOLDER/path/to/phyda.nc'
TARGET_MONTHS <- c(6, 7, 8) # Ej: JJA = 6,7,8
OUTPUT_TXT  <- 'PLACEHOLDER/path/to/phyda_annual.txt'

# Coordenadas (Javalambre / Teruel)
LAT_CENTER  <- 40.11
LON_CENTER  <- -1.04 

# Caja para promedio espacial
LAT_MIN     <- 39.5 
LAT_MAX     <- 40.7
LON_MIN_360 <- 358.0 
LON_MAX_360 <- 359.8
# ------------------------------


# ==============================================================================
# 2. FUNCIÓN DE EXTRACCIÓN MEJORADA
# ==============================================================================

process_phyda_smart <- function(nc_path, lat_range, lon_range_360, center_coords, months) {
  
  nc <- nc_open(nc_path)
  
  message(paste("Abriendo:", basename(nc_path)))
  
  # --- A. GESTIÓN DE COORDENADAS ---
  lats <- ncvar_get(nc, "lat")
  lons <- ncvar_get(nc, "lon") # Puede ser 0-360 o -180-180
  years <- ncvar_get(nc, "time")
  
  # Detectar si el NetCDF usa 0-360 o -180/180
  use_0_360 <- any(lons > 180)
  
  # Ajustar la longitud objetivo según el formato del archivo
  if (use_0_360) {
    target_lon_min <- lon_range_360[1]
    target_lon_max <- lon_range_360[2]
    center_lon_adj <- ifelse(center_coords$lon < 0, center_coords$lon + 360, center_coords$lon)
  } else {
    # Convertir tus 358 a -2 si el archivo usa formato -180/180
    target_lon_min <- ifelse(lon_range_360[1] > 180, lon_range_360[1] - 360, lon_range_360[1])
    target_lon_max <- ifelse(lon_range_360[2] > 180, lon_range_360[2] - 360, lon_range_360[2])
    center_lon_adj <- center_coords$lon
  }
  
  # --- B. ENCONTRAR ÍNDICES (Lógica Robusta) ---
  
  # 1. Intentar encontrar índices DENTRO de la caja
  lat_idxs <- which(lats >= lat_range[1] & lats <= lat_range[2])
  lon_idxs <- which(lons >= target_lon_min & lons <= target_lon_max)
  
  # 2. FAILSAFE: Si la caja es muy pequeña y no cae ningún centro de pixel dentro...
  # Buscamos el píxel más cercano al centro (Nearest Neighbor)
  if (length(lat_idxs) == 0) {
    message("AVISO: Ninguna latitud cae estrictamente en la caja. Usando la más cercana.")
    lat_idxs <- which.min(abs(lats - center_coords$lat))
  }
  
  if (length(lon_idxs) == 0) {
    message("AVISO: Ninguna longitud cae estrictamente en la caja. Usando la más cercana.")
    lon_idxs <- which.min(abs(lons - center_lon_adj))
  }
  
  message(paste("Píxeles seleccionados -> Lat:", length(lat_idxs), "| Lon:", length(lon_idxs)))
  message(paste("Rango Lat real:", round(min(lats[lat_idxs]),2), "-", round(max(lats[lat_idxs]),2)))
  
  # --- C. PREPARAR LECTURA OPTIMIZADA (start/count) ---
  # Esto evita leer todo el mapa mundial
  start_lon <- min(lon_idxs)
  count_lon <- length(lon_idxs)
  start_lat <- min(lat_idxs)
  count_lat <- length(lat_idxs)
  
  # Dataframe base
  df_out <- data.frame(Year = years)
  
  # --- D. VARIABLES ESPACIALES ---
  spatial_vars <- c("tas_mn", "pdsi_mn", "spei_mn")
  
  for (var in spatial_vars) {
    if (var %in% names(nc$var)) {
      message(paste("Procesando:", var))
      
      # Leemos SOLO el trozo que necesitamos
      # Asumimos dimensiones: [Lon, Lat, Time]
      subset_data <- ncvar_get(nc, var, 
                               start = c(start_lon, start_lat, 1), 
                               count = c(count_lon, count_lat, -1)) # -1 lee todos los tiempos
      
      # Calcular media espacial
      if (length(dim(subset_data)) == 3) {
        # Si hay múltiples píxeles (caja grande)
        regional_mean <- apply(subset_data, 3, mean, na.rm=TRUE)
      } else if (is.matrix(subset_data)) {
        # Caso raro: solo 1 lat y varias lon (o viceversa), o matriz tiempo x 1 pixel
        # Generalmente si count_lon=1 y count_lat=1, ncvar_get devuelve un vector (Time)
        regional_mean <- colMeans(subset_data, na.rm=TRUE) 
      } else {
        # Es un vector simple (1 pixel espacial, N tiempos)
        regional_mean <- subset_data
      }
      
      df_out[[paste0(var, "_local")]] <- regional_mean
    }
  }
  
  # --- E. ÍNDICES GLOBALES (Sin cambios, son vectores rápidos) ---
  indices_annual <- c("amo_mn", "gmt_mn", "Atl_mn", "Pac130_mn")
  for (var in indices_annual) {
    if (var %in% names(nc$var)) df_out[[var]] <- ncvar_get(nc, var)
  }
  
  # --- F. ÍNDICES EL NIÑO (Mensuales -> Anuales) ---
  indices_monthly <- c("Nino_3.4_mn", "Nino_12_mn", "Nino_3_mn", "Nino_4_mn")
  for (var in indices_monthly) {
    if (var %in% names(nc$var)) {
      raw_monthly <- ncvar_get(nc, var) 
      # Convertir a matriz 12 x Años
      # NOTA: Asegurarse de que raw_monthly es divisible por 12. PHYDA suele serlo.
      n_years_idx <- length(raw_monthly) / 12
      mat_monthly <- matrix(raw_monthly, nrow=12, ncol=n_years_idx, byrow=FALSE)
      
      if (length(months) > 1) {
        seasonal_mean <- colMeans(mat_monthly[months, ], na.rm=TRUE)
      } else {
        seasonal_mean <- mat_monthly[months, ]
      }
      # Recortar o rellenar si los años de los índices difieren del NetCDF principal
      # (Generalmente PHYDA cuadra, pero por seguridad tomamos los primeros N años)
      df_out[[var]] <- seasonal_mean[1:length(years)]
    }
  }
  
  nc_close(nc)
  return(df_out)
}

# ==============================================================================
# 3. EJECUTAR
# ==============================================================================

datos_clima <- process_phyda_smart(
  nc_path = NC_PATH, 
  lat_range = c(LAT_MIN, LAT_MAX), 
  lon_range_360 = c(LON_MIN_360, LON_MAX_360),
  center_coords = list(lat = LAT_CENTER, lon = LON_CENTER),
  months = TARGET_MONTHS
)

# Chequeo rápido
print(head(datos_clima))

# Guardar
write.table(datos_clima, OUTPUT_TXT, sep="\t", row.names=FALSE, quote=FALSE)
message("¡Listo! Archivo generado: ", OUTPUT_TXT)