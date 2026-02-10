# ==============================================================================
# 07_format_terra_txt.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   TerraClimate Raw Text Formatter.
#   Post-processing script for raw TerraClimate exports (often "long format").
#   It standardizes dates, coordinates, and ensures the data corresponds exactly 
#   to a target sampling site using nearest-neighbor logic if the export contains 
#   regional data.
#
#   Key Steps:
#   - Date Standardisation: Converts various date formats to YYYYMMDD.
#   - Spatial Filtering: Selects only the grid point closest to the Sample Lat/Lon.
#   - Reshaping: Converts valid data into the standard Year x Month matrix.
#
#   Inputs:
#   - Raw TerraClimate text file.
#   - Target Site Coordinates.
#
#   Outputs:
#   - Overwrites input file with formatted matrix.
# ==============================================================================

# --- PARÁMETROS A MODIFICAR ---
INPUT_FILE  <- "PLACEHOLDER/path/to/input_terra_data.txt"
LAT_MUESTRA <- 39.901
LON_MUESTRA <- 16.228
# ------------------------------

# Cargar librerías
library(data.table)
library(geosphere)




# Función para cambiar los nombres de las columnas
cambiar_nombres_columnas <- function(data) {
  setnames(data, old = names(data), new = c("date", "latitude", "longitude", "precipitation"))
}

# Función principal
procesar_datos <- function(input_file, lat_muestra, lon_muestra) {
  # Cargar los datos
  data <- fread(input_file, sep = "\t")
  
  # Cambiar los nombres de las columnas
  cambiar_nombres_columnas(data)
  
  # Crear una función para generar las fechas
  generate_dates <- function(start_year, start_month, n) {
    dates <- seq(as.Date(paste(start_year, start_month, "01", sep = "-")), by = "month", length.out = n)
    return(dates)
  }
  
  # Obtener las fechas únicas y generar las nuevas fechas
  unique_dates <- unique(data$date)
  new_dates <- generate_dates(1958, 1, length(unique_dates))
  
  # Crear un dataframe para mapear las fechas originales con las nuevas fechas
  date_mapping <- data.frame(original = unique_dates, new = new_dates)
  
  # Reemplazar las fechas en el dataframe original
  data$date <- as.Date(sapply(data$date, function(x) date_mapping$new[date_mapping$original == x]))
  
  # Redondear las columnas de latitud y longitud a dos decimales
  data$latitude <- round(data$latitude, 2)
  data$longitude <- round(data$longitude, 2)
  
  # Cambiar el formato de la fecha para eliminar los guiones
  data$date <- format(data$date, "%Y%m%d")
  
  # Filtrar las filas más cercanas geográficamente para la primera fecha
  datos_primera_fecha <- data[date == min(date)]
  datos_muestra <- datos_primera_fecha[which.min(distHaversine(cbind(longitude, latitude), c(lon_muestra, lat_muestra))), ]
  
  # Obtener las coordenadas más cercanas
  lat_cercana <- datos_muestra$latitude
  lon_cercana <- datos_muestra$longitude
  
  # Filtrar los datos usando las coordenadas más cercanas
  datos_muestra <- data[latitude == lat_cercana & longitude == lon_cercana]
  
  # Agrupar por año y mes y sumar valores de precipitation
  datos_muestra[, `:=`(year = as.numeric(substr(date, 1, 4)), month = as.numeric(substr(date, 5, 6)))]
  datos_muestra <- datos_muestra[, .(precipitation_sum = sum(precipitation, na.rm = TRUE)), by = .(year, month)]
  
  # Reestructurar en formato filas por año, columnas por mes
  datos_muestra_wide <- dcast(datos_muestra, year ~ month, value.var = "precipitation_sum")
  
  # Guardar el archivo modificado, sobrescribiendo el original sin encabezado
  fwrite(datos_muestra_wide, file = input_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Llamar a la función con el archivo de entrada y las coordenadas
procesar_datos(INPUT_FILE, LAT_MUESTRA, LON_MUESTRA)
