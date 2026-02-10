# ==============================================================================
# 05_correlacion_espacial_eobs.R
# Autor: Marcos Marín-Martín
# Fecha: 2026-02-09
# Descripción: Correlación espacial entre cronologías y precipitación acumulada (NetCDF E-OBS).
#              Genera mapas de correlación y destaca el píxel de máxima correlación.
# ==============================================================================

# --- PARÁMETROS A MODIFICAR ---
# --- PARÁMETROS A MODIFICAR ---
NC_FILE              <- "PLACEHOLDER/path/to/eobs01.nc"
CHRONOLOGY_FILE_PATH <- "PLACEHOLDER/path/to/chronology.txt"
START_YEAR           <- 1950
END_YEAR             <- 2022
ACCUMULATED_DAYS     <- 320 # Días para acumulación de precipitación
END_MONTH            <- 7   # Mes de finalización (Julio)
END_DAY              <- 1   # Día de finalización
PLOT_FILENAME_PNG    <- "correlacion_espacial_terra_final.png" 
PLOT_FILENAME_PDF    <- "correlacion_espacial_terra_final.pdf"
# ------------------------------

# --- 1. Gestión de Paquetes ---
library(ncdf4)
library(terra)
library(lubridate)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(PCICt)

# --- 1. Gestión de Paquetes ---

cat("Paquetes necesarios cargados.\n")


# --- 3. Abrir archivo NetCDF (SIN leer 'rr' aún y SIN cerrarlo) ---
tryCatch({
  # Dejamos nc_data abierto para leer dentro del bucle
  nc_data <- nc_open(NC_FILE) 
}, error = function(e) {
  stop("Error abriendo el archivo NetCDF: ", NC_FILE, "\n", e$message)
})

# Extraer dimensiones y metadatos (lon, lat, time_vals, units, calendar)
lon <- tryCatch(ncvar_get(nc_data, "longitude"), error = function(e) stop("Variable 'longitude' no encontrada."))
lat <- tryCatch(ncvar_get(nc_data, "latitude"), error = function(e) stop("Variable 'latitude' no encontrada."))
time_vals <- tryCatch(ncvar_get(nc_data, "time"), error = function(e) stop("Variable 'time' no encontrada."))
time_units <- tryCatch(ncatt_get(nc_data, "time", "units")$value, error = function(e) "Error obteniendo unidades de tiempo.")
calendar <- tryCatch(ncatt_get(nc_data, "time", "calendar")$value, error = function(e) {warning("Atributo 'calendar' no encontrado, asumiendo gregoriano."); "gregorian"})

# Obtener dimensiones directamente de la info del NetCDF para robustez
n_lon <- length(lon)
n_lat <- length(lat)
# NO LEEMOS rr_array aquí para ahorrar memoria

cat("Metadatos NetCDF leídos.\n")

# --- 4. Convertir tiempo a fechas ---
cat("Convirtiendo fechas...\n") # Mensaje de inicio
tryCatch({ 
  # --- Todo el código que estaba dentro del 'try' original va aquí ---
  origin_parts <- strsplit(time_units, " since ")[[1]]
  time_unit_word <- tolower(strsplit(origin_parts[1], " ")[[1]][1]) 
  origin_timestamp_str <- origin_parts[2]
  
  if (!is.null(calendar) && !(tolower(calendar) %in% c("gregorian", "standard", "proleptic_gregorian"))) {
    cat("Usando calendario no estándar:", calendar, " con PCICt.\n")
    origin_pcict <- as.PCICt(origin_timestamp_str, cal = calendar)
    # Convertir unidades a segundos para PCICt
    if (time_unit_word == "days") {
      dates_pcict <- origin_pcict + time_vals * 86400
    } else if (time_unit_word == "hours") {
      dates_pcict <- origin_pcict + time_vals * 3600
    } else if (time_unit_word == "seconds") {
      dates_pcict <- origin_pcict + time_vals
    } else {
      stop("Unidades de tiempo no soportadas para PCICt: ", time_unit_word)
    }
    dates <- as.Date(dates_pcict) # Convertir a objeto Date para comparación fácil
  } else {
    cat("Usando calendario estándar con lubridate.\n")
    origin_date <- as.POSIXct(origin_timestamp_str, tz = "UTC")
    if (time_unit_word == "days") {
      dates <- as.Date(origin_date + lubridate::days(time_vals))
    } else if (time_unit_word == "hours") {
      dates <- as.Date(origin_date + lubridate::hours(time_vals))
    } else if (time_unit_word == "seconds") {
      dates <- as.Date(origin_date + lubridate::seconds(time_vals))
    } else {
      stop("Unidades de tiempo no soportadas para lubridate: ", time_unit_word)
    }
  }
  # --- Fin del código del 'try' original ---
  
  cat("Fechas convertidas correctamente.\n") # Mensaje si todo fue bien
  
}, error = function(e) { 
  # --- El código que estaba en el 'catch' original va aquí, dentro de la función 'error' ---
  stop("Error convirtiendo las fechas: ", e$message) 
}) # Cierre del tryCatch
# Ya no necesitamos el cat("Fechas convertidas.\n") aquí fuera, se maneja dentro o se detiene.

# --- 5. Calcular la precipitación acumulada (Leyendo 'rr' por año) ---
years <- seq(START_YEAR, END_YEAR)
n_years <- length(years)

# Inicializar array para resultados [lon, lat, year]
accumulated_rr_array <- array(NA, dim = c(n_lon, n_lat, n_years),
                              dimnames = list(lon = lon, lat = lat, year = years))

cat("Calculando precipitación acumulada anual (leyendo por año)...\n")
pb <- txtProgressBar(min = 0, max = n_years, style = 3) # Barra de progreso
for (i in 1:n_years) {
  current_year <- years[i]
  
  # Definir fechas de inicio y fin
  end_date <- as.Date(sprintf("%d-%02d-%02d", current_year, END_MONTH, END_DAY))
  start_date <- end_date - lubridate::days(ACCUMULATED_DAYS) 
  
  # Encontrar índices de tiempo para la ventana del año actual
  mask_indices <- which(dates >= start_date & dates < end_date)
  
  if (length(mask_indices) > 0) {
    # --- Lectura de 'rr' solo para los índices de este año ---
    # Encontrar el índice mínimo y máximo para leer un bloque contiguo
    min_idx <- min(mask_indices)
    max_idx <- max(mask_indices)
    
    # Definir el inicio (start) y la cuenta (count) para ncvar_get
    # Formato: [lon, lat, time]
    start_vec <- c(1, 1, min_idx) 
    # Leer todos los lon, todos los lat, y el bloque de tiempo desde min_idx hasta max_idx
    count_vec <- c(n_lon, n_lat, max_idx - min_idx + 1) 
    
    # Leer el CHUNK de datos que contiene los días necesarios para este año
    rr_chunk_year <- tryCatch({
      ncvar_get(nc_data, "rr", start = start_vec, count = count_vec)
    }, error = function(e){
      # Advertir si falla la lectura para este año y continuar
      warning(sprintf("Error leyendo chunk NetCDF para año %d (índices %d a %d): %s. Saltando año.", 
                      current_year, min_idx, max_idx, e$message))
      NULL # Devolver NULL para indicar fallo
    })
    
    # Si la lectura fue exitosa (no es NULL)
    if (!is.null(rr_chunk_year)){
      # Calcular los índices RELATIVOS dentro del chunk que corresponden a mask_indices
      # Ej: si mask_indices = [10, 11, 15] y min_idx = 10, los índices en el chunk son [1, 2, 6]
      indices_in_chunk <- mask_indices - min_idx + 1
      
      # Subseleccionar SÓLO los días exactos requeridos del chunk leído
      # Usar drop = FALSE por si indices_in_chunk tiene longitud 1
      rr_subset <- rr_chunk_year[, , indices_in_chunk, drop = FALSE] 
      
      # Sumar sobre la dimensión de tiempo (la tercera dimensión) del subset
      accumulated_rr_array[, , i] <- apply(rr_subset, MARGIN = c(1, 2), FUN = sum, na.rm = TRUE)
      
      # Opcional: Liberar memoria explícitamente (puede o no ser necesario/útil)
      # rm(rr_chunk_year, rr_subset)
      # gc() # Llamar al garbage collector
      
    } else {
      # Si falló la lectura del chunk, asignar NA a todo el año
      accumulated_rr_array[, , i] <- NA 
    }
    # --- Fin de la lectura de 'rr' para este año ---
    
  } else {
    # Si no hay datos en el período para este año, llenar con NA 
    accumulated_rr_array[, , i] <- NA 
  }
  setTxtProgressBar(pb, i) # Actualizar barra de progreso
}
close(pb)

# --- IMPORTANTE: Cerrar el archivo NetCDF AHORA que terminamos de leer ---
nc_close(nc_data) 
cat("\nCálculo de precipitación acumulada completado. Archivo NetCDF cerrado.\n")

# --- 6. Cargar archivo de cronología ---
cat("Cargando archivo de cronología...\n")
tryCatch({
  # --- Corrección aquí ---
  chronology_data <- read.table(
    CHRONOLOGY_FILE_PATH, 
    header = TRUE,        # CORRECTO: La primera línea es el encabezado
    # skip = 1,           # INCORRECTO: No saltar el encabezado
    # header = FALSE,     # INCORRECTO: Sí hay encabezado
    na.strings = "NA",    # CORRECTO: Manejar "NA" como NA de R
    sep = "\t",           # Especificar tabulador como separador (más seguro)
    # col.names = c("Year", "Var2", "Residuals") # INCORRECTO: No necesario y erróneo
    check.names = FALSE   # Opcional: previene que R cambie nombres como 'samp.depth' a 'samp..depth'
  ) 
  # --- Fin Corrección ---
  
  # Verificar que se leyó algo
  if (nrow(chronology_data) == 0) {
    stop("El archivo de cronología se leyó pero está vacío.")
  }
  # Verificar que las columnas esperadas existen (con los nombres del archivo)
  if (!("year" %in% names(chronology_data)) || !("res" %in% names(chronology_data))) {
    stop("Las columnas 'year' o 'res' no se encontraron en el archivo de cronología. Nombres encontrados: ", 
         paste(names(chronology_data), collapse=", "))
  }
  
  
}, error = function(e) {
  stop("Error leyendo el archivo de cronología: ", CHRONOLOGY_FILE_PATH, "\n", e$message)
})

# Filtrar los datos de cronología para los años de interés
# --- Corrección aquí: Usar 'year' (minúscula) ---
chronology_filtered <- subset(chronology_data, year >= START_YEAR & year <= END_YEAR)

# Asegurar que la cronología cubra todos los años y esté ordenada
year_df <- data.frame(Year = years) # Mantener 'Year' aquí porque 'years' es la secuencia
# --- Corrección aquí: Usar 'year' (minúscula) en 'by.x' para el merge ---
chronology_aligned <- merge(year_df, chronology_filtered, 
                            by.x = "Year", by.y = "year",  # Especificar columnas de unión
                            all.x = TRUE)

# Extraer los residuos alineados
# --- Corrección aquí: Usar 'res' (minúscula) ---
filtered_chronology_residuals <- chronology_aligned$res

# Validar longitud (después del merge, deberían coincidir)
if (length(filtered_chronology_residuals) != n_years) {
  stop("Error inesperado: La longitud de la cronología alineada no coincide con el número de años.")
}
if(sum(!is.na(filtered_chronology_residuals)) < 3) {
  stop("La cronología filtrada para el período ", START_YEAR, "-", END_YEAR, 
       " tiene menos de 3 valores no-NA en la columna 'res'. No se puede calcular la correlación.")
}
cat("Datos de cronología cargados, filtrados y alineados correctamente.\n")


# --- 7. Calcular correlaciones espaciales y R ---
# (Esta sección permanece igual)
correlations <- matrix(NA, nrow = n_lon, ncol = n_lat)
p_values <- matrix(NA, nrow = n_lon, ncol = n_lat)
cat("Calculando correlaciones espaciales (puede tardar)... \n")
pb <- txtProgressBar(min = 0, max = n_lon * n_lat, style = 3)
progress_count <- 0
for (i in 1:n_lon) {
  for (j in 1:n_lat) {
    pixel_ts <- accumulated_rr_array[i, j, ]
    valid_indices <- !is.na(filtered_chronology_residuals) & !is.na(pixel_ts)
    n_valid <- sum(valid_indices)
    if (n_valid > 2) { 
      test_result <- tryCatch({
        chrono_valid <- filtered_chronology_residuals[valid_indices]
        pixel_valid <- pixel_ts[valid_indices]
        if(sd(chrono_valid, na.rm=TRUE) == 0 | sd(pixel_valid, na.rm=TRUE) == 0) {
          list(estimate = NA, p.value = NA) 
        } else { cor.test(chrono_valid, pixel_valid, method = "pearson") }
      }, error = function(e) { list(estimate = NA, p.value = NA) })
      correlations[i, j] <- test_result$estimate
      p_values[i, j] <- test_result$p.value
    } 
    progress_count <- progress_count + 1
    setTxtProgressBar(pb, progress_count)
  }
}
close(pb)
cat("\nCálculo de correlaciones completado.\n")

# --- 8. Filtrar R significativos (p < 0.05) y positivos ---
# (Esta sección permanece igual)
significant_r <- correlations 
significant_r[is.na(p_values) | p_values >= 0.05 | is.na(correlations) | correlations <= 0] <- NA
cat(sprintf("Número de píxeles con R significativo y positivo: %d\n", sum(!is.na(significant_r))))

# --- 9. Crear la paleta de colores ---
# (Esta sección permanece igual)
colors_hex <- c('#F1BC9C', '#E18466', '#D24734', '#C40A0C')
custom_cmap <- colorRampPalette(colors_hex)

# --- 10. Plotear el mapa usando 'terra' y EXPORTAR ---

# Determinar si lat es ascendente o descendente
lat_order_asc <- all(diff(lat) > 0)
if (lat_order_asc) { cat("Orden de latitud original detectado: Ascendente (Sur a Norte)\n") 
} else { cat("Orden de latitud original detectado: Descendente (Norte a Sur)\n") }

# Crear el SpatRaster (lógica de inversión ya corregida)
tryCatch({
  r_transposed <- t(significant_r) 
  if (lat_order_asc) {
    cat("Invirtiendo filas de la matriz transpuesta para que coincidan con la expectativa de terra (N->S).\n")
    r_transposed <- r_transposed[nrow(r_transposed):1, , drop = FALSE] 
  } else {
    cat("La matriz transpuesta ya tiene filas N->S. No se necesita inversión para terra.\n")
  }
  terra_r <- terra::rast(r_transposed, crs = "EPSG:4326") 
  cat(sprintf("Estableciendo extensión: xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f\n", 
              min(lon), max(lon), min(lat), max(lat)))
  terra::ext(terra_r) <- c(min(lon), max(lon), min(lat), max(lat))
}, error = function(e){
  stop("Error creando el SpatRaster: ", e$message)
})


# Obtener datos de países/costas (sin cambios)
world_sf <- ne_countries(scale = "medium", returnclass = "sf")
world_sf <- sf::st_make_valid(world_sf) 
world_sf <- sf::st_transform(world_sf, crs = "EPSG:4326")
world_vect <- terra::vect(world_sf) 

# --- Crear una FUNCIÓN con todos los comandos de ploteo ---
# Esto evita repetir el código para PNG y PDF
plot_correlation_map <- function() {
  # Plot base del raster
  terra::plot(terra_r, 
              col = custom_cmap(100), 
              type = "continuous",    
              main = "Correlación espacial con precipitación acumulada (R significativa > 0, p < 0.05)",
              axes = FALSE, 
              mar = c(4, 4, 3, 1), 
              plg = list(title = "R", title.cex=0.9, cex=0.8) 
  )
  
  # Añadir costas y fronteras
  terra::plot(world_vect, add = TRUE, border = "gray50", lwd = 0.6) 
  
  # Añadir Escala Gráfica
  terra::sbar(d = 500, xy = "bottomleft", type = "bar", divs = 2, 
              below = "km", cex = 0.8, lwd = 1.5, lonlat = TRUE)  
  
  # Añadir Flecha Norte
  usr <- par("usr") 
  arrow_x <- usr[1] + (usr[2] - usr[1]) * 0.95 
  arrow_y_base <- usr[3] + (usr[4] - usr[3]) * 0.90 
  arrow_y_tip <- usr[3] + (usr[4] - usr[3]) * 0.95 
  arrows(x0 = arrow_x, y0 = arrow_y_base, x1 = arrow_x, y1 = arrow_y_tip, 
         length = 0.08, angle = 30, lwd = 1.5)     
  text(x = arrow_x, y = arrow_y_tip, label = "N", pos = 3, offset = 0.5, cex = 0.9)       
}
# --- Fin de la función de ploteo ---


# --- Exportar a PNG ---
cat(sprintf("Exportando a PNG: %s\n", PLOT_FILENAME_PNG))
png(PLOT_FILENAME_PNG, width = 1200, height = 900, res = 100) 
# Llamar a la función de ploteo
plot_correlation_map() 
dev.off() # Cerrar dispositivo PNG
cat(sprintf("Mapa PNG guardado como: %s\n", PLOT_FILENAME_PNG))


# --- Exportar a PDF Vectorial ---
cat(sprintf("Exportando a PDF: %s\n", PLOT_FILENAME_PDF))
# Para PDF, width y height son en pulgadas. Ajusta según necesidad.
# 1200px / 100 dpi = 12 pulgadas; 900px / 100 dpi = 9 pulgadas
pdf(PLOT_FILENAME_PDF, width = 10, height = 7.5, useDingbats=FALSE) # useDingbats=FALSE es bueno para compatibilidad
# Llamar a la MISMA función de ploteo
plot_correlation_map() 
dev.off() # Cerrar dispositivo PDF
cat(sprintf("Mapa PDF guardado como: %s\n", PLOT_FILENAME_PDF))



# --- 11. Encontrar el valor más alto de R y sus coordenadas ---
# (Esta sección permanece igual)
max_r_value <- max(significant_r, na.rm = TRUE)
if (is.finite(max_r_value)) {
  max_indices_matrix <- which(significant_r == max_r_value, arr.ind = TRUE)
  max_lon_index <- max_indices_matrix[1, 1]
  max_lat_index <- max_indices_matrix[1, 2]
  max_r_longitude <- lon[max_lon_index]
  max_r_latitude <- lat[max_lat_index]
  cat(sprintf("El valor más alto de R significativo y positivo es %.4f\n", max_r_value))
  cat(sprintf("Se encuentra en las coordenadas (longitud, latitud): (%.4f, %.4f)\n", max_r_longitude, max_r_latitude))
} else {
  cat("No se encontraron valores de R significativos y positivos finitos.\n")
}

cat("Análisis completado.\n")