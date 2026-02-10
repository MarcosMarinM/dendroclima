# ==============================================================================
# 06_format_rocio_txt.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Precipitation Data Processor & Aggregator (Daily to Monthly).
#   Converts raw daily precipitation data (typically from local station logs or 
#   specific gridded text exports) into clean, analytical formats.
#
#   Capabilities:
#   1. Daily Aggregation: Averages precipitation measurements for the same day 
#      (handling multiple sensors/entries).
#   2. Monthly Upscaling: Sums daily values to compute Total Monthly Precipitation.
#   3. Matrix Generation: Reshapes monthly data into a Wide Format (Year x Month) 
#      suitable for correlation analysis.
#
#   Inputs:
#   - Raw TSV/TXT file with daily records.
#
#   Outputs:
#   - Clean Daily Average file.
#   - Monthly Total Matrix file.
# ==============================================================================

# --- PARÁMETROS A MODIFICAR ---
INPUT_FILE          <- "PLACEHOLDER/path/to/input_precip.txt"
OUTPUT_FILE_DAILY   <- "PLACEHOLDER/path/to/output_daily.txt"
OUTPUT_FILE_MONTHLY <- "PLACEHOLDER/path/to/output_monthly.txt"

START_DATE_STR      <- "1951-01-01"

# Mapeo de columnas
INPUT_COL_NAMES_ORIGINAL <- c("time", "height", "rlat", "rlon", "precipitation")
INPUT_COL_NAMES_NEW      <- c("time", "height", "latitude", "longitude", "precipitation")
# ------------------------------


# Verifica que los nombres originales y nuevos tengan la misma longitud
if(length(INPUT_COL_NAMES_ORIGINAL) != length(INPUT_COL_NAMES_NEW)){
  stop("Error de configuración: 'INPUT_COL_NAMES_ORIGINAL' e 'INPUT_COL_NAMES_NEW' deben tener el mismo número de elementos.")
}
#-------------------------------------------------------------------------------#

# --- INICIO DEL CÓDIGO ---

# Cargar librerías
library(data.table)

# Función principal para procesar datos y generar ambos outputs
procesar_datos_completo <- function(input_path, output_path_daily, output_path_monthly, start_date, original_colnames, new_colnames) {
  
  # Validar que el archivo de entrada existe
  if (!file.exists(input_path)) {
    stop("Error Crítico: El archivo de entrada no existe en la ruta especificada: ", input_path)
  }
  
  message("Iniciando procesamiento del archivo: ", input_path)
  message("-> Se generará un archivo de salida diario en: ", output_path_daily)
  message("-> Se generará un archivo de salida mensual en: ", output_path_monthly)
  message("-> El archivo de entrada original NO será modificado.")
  
  # Cargar los datos (asumiendo encabezado, como en el ejemplo y script original)
  tryCatch({
    data <- fread(input_path, sep = "\t", header = TRUE, check.names = FALSE)
    message("Archivo de entrada cargado correctamente.")
    # Validar que las columnas esperadas (originales) existan
    if (!all(original_colnames %in% names(data))) {
      missing_cols <- original_colnames[!original_colnames %in% names(data)]
      stop("Error: Faltan las siguientes columnas esperadas en el archivo de entrada: ", paste(missing_cols, collapse=", "))
    }
  }, error = function(e) {
    stop("Error al cargar el archivo de entrada '", input_path, "': ", conditionMessage(e))
  })
  
  # Renombrar columnas (como en el script original)
  setnames(data, old = original_colnames, new = new_colnames)
  message("Columnas renombradas a: ", paste(new_colnames, collapse=", "))
  
  # Validar columnas clave DESPUÉS de renombrar
  if (!"time" %in% names(data) || !"precipitation" %in% names(data)) {
    stop("Error: Las columnas 'time' o 'precipitation' no se encontraron después de renombrar.")
  }
  
  # Reemplazar NA/NaN en precipitación por 0 (como en el script original)
  if (!is.numeric(data$precipitation)) {
    message("Columna 'precipitation' no es numérica. Intentando convertir...")
    data[, precipitation := suppressWarnings(as.numeric(gsub(",", ".", as.character(precipitation))))]
    if(any(is.na(data$precipitation))) {
      warning("Algunos valores de precipitación no pudieron ser convertidos a número y se tratarán como NA.")
    }
  }
  num_nas_before = sum(is.na(data$precipitation))
  data[is.na(precipitation), precipitation := 0]
  message(paste("Valores NA/NaN en 'precipitation' reemplazados por 0. Total reemplazados:", num_nas_before - sum(is.na(data$precipitation))))
  
  
  # Generar fechas (como en el script original)
  generate_dates <- function(start_dt, n) {
    start_dt_obj <- tryCatch(as.Date(start_dt), error = function(e) stop("Error: 'start_date_str' ('", start_dt,"') no es válida."))
    start_dt_obj + seq(0, by = 1, length.out = n)
  }
  
  unique_times <- unique(data$time) # Usar unique() como en el original
  n_unique_times <- length(unique_times)
  if (n_unique_times == 0) stop("Error: No se encontraron valores únicos en 'time'.")
  message(paste("Número de valores de tiempo ('time') únicos encontrados:", n_unique_times))
  fechas <- generate_dates(start_date, n_unique_times)
  date_mapping <- data.table(time = unique_times, date = fechas)
  
  # Mapear fechas (como en el script original)
  data <- merge(data, date_mapping, by = "time", all.x = TRUE, sort = FALSE)
  message("Fechas mapeadas a los datos.")
  if (any(is.na(data$date))) {
    warning(sum(is.na(data$date)), " registros no pudieron ser mapeados a una fecha y serán excluidos.")
    data <- data[!is.na(date)]
    if (nrow(data) == 0) stop("Error fatal: Ningún registro mapeado a fecha.")
  }
  
  # --- PASO INTERMEDIO: Calcular la media diaria (base para ambos) ---
  message("Calculando la media diaria de precipitación...")
  # Usar `daily_mean` como nombre de columna, igual que en script original para cálculo mensual
  daily_avg <- data[, .(daily_mean = mean(precipitation, na.rm = TRUE)), keyby = date]
  message("Cálculo de media diaria completado.")
  
  # --- TAREA 1: GENERAR Y GUARDAR EL OUTPUT DIARIO (en archivo nuevo) ---
  message("Generando archivo de salida DIARIO...")
  output_data_daily <- daily_avg[, .(
    date_formatted = format(date, "%Y%m%d"),
    precipitation = round(daily_mean, 2)
  )]
  setnames(output_data_daily, "date_formatted", "date") # Renombrar para el output
  
  header_line_daily <- '"date"\t"precipitation"'
  tryCatch({
    # Escribir encabezado
    writeLines(header_line_daily, output_path_daily)
    # Añadir datos
    fwrite(output_data_daily, file = output_path_daily, sep = "\t",
           row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
    message(paste("-> Archivo de salida DIARIO generado con éxito en:", output_path_daily))
  }, error = function(e){
    warning("Error al escribir el archivo de salida diario: ", conditionMessage(e))
    # Continuar para intentar generar el mensual
  })
  
  
  # --- TAREA 2: GENERAR Y GUARDAR EL OUTPUT MENSUAL (en archivo nuevo, sin encabezado) ---
  # Replicar lógica del script original para el cálculo mensual
  message("Generando archivo de salida MENSUAL...")
  tryCatch({
    daily_avg_monthly <- copy(daily_avg) # Usar la tabla de medias diarias
    daily_avg_monthly[, `:=`(year = year(date), month = month(date))] # Añadir año/mes
    
    # Calcular suma mensual (como en script original)
    monthly_sum <- daily_avg_monthly[, .(monthly_precipitation = sum(daily_mean, na.rm = TRUE)), by = .(year, month)]
    
    # Reestructurar a formato ancho (como en script original)
    setkey(monthly_sum, year, month) # Asegurar orden para dcast
    monthly_wide <- dcast(monthly_sum, year ~ month, value.var = "monthly_precipitation")
    
    # Guardar el archivo mensual en la NUEVA ruta, SIN encabezado (como pedía el original para su output)
    fwrite(monthly_wide, file = output_path_monthly, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, na = "NA")
    message(paste("-> Archivo de salida MENSUAL generado con éxito en:", output_path_monthly))
    
  }, error = function(e) {
    warning("Error durante la generación o guardado del archivo mensual: ", conditionMessage(e))
  })
  
  message("\nProcesamiento completado.")
}

# --- FIN DEL CÓDIGO ---

# Llamar a la función principal con manejo de errores
tryCatch({
  procesar_datos_completo(INPUT_FILE, OUTPUT_FILE_DAILY, OUTPUT_FILE_MONTHLY, START_DATE_STR, INPUT_COL_NAMES_ORIGINAL, INPUT_COL_NAMES_NEW)
}, error = function(e) {
  message("\n-----------------------------------")
  message("ERROR DURANTE LA EJECUCIÓN:")
  message(conditionMessage(e)) # Mensaje de error
  message("-----------------------------------")
  # traceback() # Descomentar si necesitas ver la traza de llamadas para depurar
})