# ==============================================================================
# 02_main_detrending_analysis.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Master Script for Dendroclimatic Series Standardisation and Analysis.
#   This comprehensive pipeline handles the core processing of tree-ring width (TRW) 
#   data into standardised chronologies.
#
#   Capabilities:
#   1. Selective Detrending: Applies multiple standardisation methods (Splines of 
#      various stiffness, Modified Negative Exponential, Mean, Autoregressive, etc.) 
#      based on user configuration.
#   2. Chronology Building: Generates both Standard (STD) and Residual (RES) 
#      chronologies from the detrended series.
#   3. RWI Export: Saves individual detrended series (Ring Width Indices) in 
#      Tucson format for further analysis.
#   4. Visualisation: Produces comparative PDF plots of the different detrending 
#      outcomes to aid in method selection.
#   5. Climate Correlations (Optional): Can perform immediate correlation analysis 
#      against a climate dataset to validate the sensitivity of the generated chronologies.
#
#   Inputs:
#   - Raw RWL file (Tucson format).
#   - Climate data file (for optional analysis).
#
#   Outputs:
#   - Detrended RWI files (.rwi).
#   - Chronology summaries (.txt).
#   - Diagnostic Plots (.pdf).
#   - Chronology summaries (.txt).
#   - Diagnostic Plots (.pdf).
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. CARGA DE LIBRERÍAS
# ------------------------------------------------------------------------------
library(dplR)
library(dendroTools)
library(ggplot2)
library(dplyr)
library(reshape2)
library(colourvalues)

# --- PARÁMETROS A MODIFICAR ---

# 1. CONTROL DE EJECUCIÓN
EJECUTAR_DETRENDING    <- TRUE    # Generar cronologías y RWIs
EJECUTAR_CORRELACIONES <- FALSE   # Calcular correlaciones climáticas

# 2. RUTAS DE ARCHIVOS (Inputs)
INPUT_RWL_FILE   <- 'PLACEHOLDER/path/to/input_data.rwl'
INPUT_CLIMA_FILE <- 'PLACEHOLDER/path/to/climate_data.txt'
DIR_SALIDA_BASE  <- 'PLACEHOLDER/path/to/output_dir'

# 3. FILTRO DE AÑOS (Solo para análisis climático)
FILTRAR_ANIOS <- FALSE
ANIO_INICIO   <- 1931   
ANIO_FIN      <- 2000    

# 4. CONFIGURACIÓN DE ANÁLISIS
ALPHA_LEVEL      <- 0.05
USAR_METODOS <- list(
  Spline_32     = FALSE,
  Spline_67     = FALSE,
  Spline_100    = FALSE,
  Spline_300    = TRUE,  # Método activo por defecto
  ModNegExp     = FALSE,
  Mean          = FALSE,  
  Ar            = FALSE,
  Friedman      = FALSE,  
  ModHugershoff = FALSE,
  AgeDepSpline  = FALSE 
)
# ------------------------------


# ------------------------------------------------------------------------------
# 2. PREPARACIÓN
# ------------------------------------------------------------------------------
# (Libraries loaded at top)

# Rutas derivadas
dir_cronos_output <- file.path(DIR_SALIDA_BASE, "RWL", "Cronologias")
dir_rwi_output    <- file.path(DIR_SALIDA_BASE, "RWL", "detrended_rwi") 
dir_corr_output   <- file.path(DIR_SALIDA_BASE, "Correlaciones")

# Crear directorios
if (!dir.exists(dir_cronos_output)) dir.create(dir_cronos_output, recursive = TRUE)
if (!dir.exists(dir_rwi_output)) dir.create(dir_rwi_output, recursive = TRUE)
if (!dir.exists(dir_corr_output)) dir.create(dir_corr_output, recursive = TRUE)

nombre_rwl_base <- tools::file_path_sans_ext(basename(INPUT_RWL_FILE))


# ==============================================================================
# PARTE 1: MULTI-DETRENDING SELECTIVO (+ RWI)
# ==============================================================================
if (EJECUTAR_DETRENDING) {
  message("\n>>> INICIANDO PARTE 1: MULTI-DETRENDING SELECTIVO")
  message("    Archivo entrada: ", basename(INPUT_RWL_FILE))
  
  data <- read.rwl(INPUT_RWL_FILE)
  
  message("    Aplicando transformación de potencia (powt)...")
  data_pt <- powt(data)
  
  years <- as.numeric(rownames(data))
  x_range <- range(years, na.rm = TRUE)
  
  # Lista vacía donde iremos metiendo SOLO lo calculado
  methods_rwi <- list()
  
  message("    Calculando métodos seleccionados (TRUE)...")
  
  # --- CÁLCULO DE MÉTODOS ---
  if (isTRUE(USAR_METODOS$Spline_32))  methods_rwi$Spline_32  <- detrend(data_pt, method = "Spline", nyrs = 32)
  
  # --- AÑADIDO: Spline 67% ---
  # Al no especificar 'nyrs', dplR usa por defecto 0.67 * longitud de la serie
  if (isTRUE(USAR_METODOS$Spline_67))  methods_rwi$Spline_67  <- detrend(data_pt, method = "Spline") 
  if (isTRUE(USAR_METODOS$Spline_100)) methods_rwi$Spline_100 <- detrend(data_pt, method = "Spline", nyrs = 100)
  if (isTRUE(USAR_METODOS$Spline_300)) methods_rwi$Spline_300 <- detrend(data_pt, method = "Spline", nyrs = 300)
  if (isTRUE(USAR_METODOS$ModNegExp))  methods_rwi$ModNegExp  <- detrend(data_pt, method = "ModNegExp")
  if (isTRUE(USAR_METODOS$Mean))       methods_rwi$Mean       <- detrend(data_pt, method = "Mean")
  if (isTRUE(USAR_METODOS$Ar))         methods_rwi$Ar         <- detrend(data_pt, method = "Ar")
  if (isTRUE(USAR_METODOS$Friedman))   methods_rwi$Friedman   <- detrend(data_pt, method = "Friedman")
  if (isTRUE(USAR_METODOS$ModHugershoff)) methods_rwi$ModHugershoff <- detrend(data_pt, method = "ModHugershoff")
  if (isTRUE(USAR_METODOS$AgeDepSpline))  methods_rwi$AgeDepSpline  <- detrend(data_pt, method = "AgeDepSpline")
  
  if (length(methods_rwi) == 0) stop("ERROR: No has seleccionado ningún método en USAR_METODOS.")
  
  # Generar Cronologías (TXT) y Archivos Individuales (RWI)
  plot_data_list <- list()
  
  for (name in names(methods_rwi)) {
    rwi_matrix <- methods_rwi[[name]]
    
    # 1. Generar Cronología Maestra
    crono <- chron(rwi_matrix, prewhiten = TRUE)
    
    col_std  <- if("std" %in% names(crono)) crono$std else rep(NA, nrow(crono))
    col_res  <- if("res" %in% names(crono)) crono$res else rep(NA, nrow(crono))
    col_samp <- crono$samp.depth
    
    df_export <- data.frame(year = as.numeric(row.names(crono)), std = col_std, res = col_res, samp.depth = col_samp)
    
    # Guardar Cronología TXT
    nombre_txt <- paste0(nombre_rwl_base, "_", name, ".txt")
    write.table(df_export, file = file.path(dir_cronos_output, nombre_txt), sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Guardar RWI Individual (STD)
    nombre_rwi_std <- paste0(nombre_rwl_base, "_", name, "_STD.rwi")
    write.rwl(rwi_matrix, fname = file.path(dir_rwi_output, nombre_rwi_std), format = "tucson")
    
    # Guardar RWI Individual (RES) si es posible
    rwi_res_matrix <- rwi_matrix
    calculate_res <- TRUE
    
    tryCatch({
      for (i in 1:ncol(rwi_matrix)) {
        series <- na.omit(rwi_matrix[, i])
        if (length(series) > 5) {
          ar_fit <- ar(series, aic = TRUE, order.max = 5)
          # Reconstruir matriz residual preservando NAs
          rwi_res_matrix[!is.na(rwi_matrix[, i]), i] <- ar_fit$resid + mean(series)
        }
      }
    }, error = function(e) { calculate_res <- FALSE })
    
    if (calculate_res) {
      nombre_rwi_res <- paste0(nombre_rwl_base, "_", name, "_RES.rwi")
      write.rwl(rwi_res_matrix, fname = file.path(dir_rwi_output, nombre_rwi_res), format = "tucson")
    }
    
    plot_data_list[[name]] <- df_export
  }
  message("    Archivos TXT (Cronologías) y RWI (Individuales) guardados.")
  
  # --- PLOTS DINÁMICOS ---
  raw_series <- rowMeans(data, na.rm = TRUE)
  colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "magenta", "darkgreen")
  
  # PDF STD
  pdf(file = file.path(dir_cronos_output, paste0(nombre_rwl_base, "_COMPARATIVA_STD.pdf")), width = 8.27, height = 11.69)
  par(mfrow = c(5, 2), mar = c(4, 4, 2, 1))
  
  i <- 1
  for (name in names(plot_data_list)) {
    df <- plot_data_list[[name]]
    col_plot <- colors[(i - 1) %% length(colors) + 1] 
    plot(df$year, df$std, type = "l", col = col_plot, lwd = 0.5,
         main = paste(name, "(Std)"), xlab = "Año", ylab = "RWI", xlim = x_range)
    abline(h = 1, col = "gray50", lty = 2)
    i <- i + 1
  }
  plot(years, raw_series, type = "l", col = "black", lwd = 0.5,
       main = "RAW Series (mm)", xlab = "Año", ylab = "mm", xlim = x_range)
  lines(lowess(years, raw_series, f = 0.1), col = "red", lwd = 0.5)
  dev.off()
  
  # PDF RES
  pdf(file = file.path(dir_cronos_output, paste0(nombre_rwl_base, "_COMPARATIVA_RES.pdf")), width = 8.27, height = 11.69)
  par(mfrow = c(5, 2), mar = c(4, 4, 2, 1))
  
  i <- 1
  for (name in names(plot_data_list)) {
    df <- plot_data_list[[name]]
    col_plot <- colors[(i - 1) %% length(colors) + 1]
    if (all(is.na(df$res))) {
      plot(1, type="n", main=paste(name, "(Res) - Sin datos"), axes=FALSE)
      text(1, 1, "No se pudo calcular AR")
    } else {
      plot(df$year, df$res, type = "l", col = col_plot, lwd = 0.5,
           main = paste(name, "(Res)"), xlab = "Año", ylab = "RWI (Res)", xlim = x_range)
      abline(h = 1, col = "gray50", lty = 2)
    }
    i <- i + 1
  }
  plot(years, raw_series, type = "l", col = "black", lwd = 0.5,
       main = "RAW Series (Ref)", xlab = "Año", ylab = "mm", xlim = x_range)
  lines(lowess(years, raw_series, f = 0.1), col = "red", lwd = 0.5)
  dev.off()
  
  message(">>> PARTE 1 COMPLETADA.\n")
}


# ==============================================================================
# PARTE 2: CORRELACIONES CLIMÁTICAS MASIVAS (FILTRADAS)
# ==============================================================================
if (EJECUTAR_CORRELACIONES) {
  
  message("\n>>> INICIANDO PARTE 2: CORRELACIONES CLIMÁTICAS")
  
  # 1. PROCESAMIENTO CLIMA
  message("    Procesando datos climáticos...")
  datos_clima <- read.table(INPUT_CLIMA_FILE, header = TRUE, sep = "\t")
  colnames(datos_clima) <- gsub('"', '', colnames(datos_clima))
  
  cols_disponibles <- names(datos_clima)
  col_variable <- cols_disponibles[!cols_disponibles %in% c("date")][1]
  message("    Variable detectada: ", col_variable)
  datos_clima$valor_analisis <- datos_clima[[col_variable]]
  
  # --- DETECCIÓN PERIODICIDAD ---
  sample_date_str <- as.character(datos_clima$date[1])
  n_char <- nchar(sample_date_str)
  TIPO_DATOS <- "DESCONOCIDO"
  
  if (n_char == 6) {
    TIPO_DATOS <- "MENSUAL"
  } else if (n_char == 8) {
    check_n <- min(nrow(datos_clima), 20)
    fechas_check <- as.Date(as.character(datos_clima$date[1:check_n]), format="%Y%m%d")
    diff_days <- median(as.numeric(diff(fechas_check)), na.rm = TRUE)
    if (diff_days > 25) TIPO_DATOS <- "MENSUAL" else TIPO_DATOS <- "DIARIO"
  } else {
    stop("ERROR: Formato de fecha desconocido.")
  }
  
  message("    >>> MODO DETECTADO: ", TIPO_DATOS)
  
  # --- CONSTRUCCIÓN MATRIZ ---
  if (TIPO_DATOS == "MENSUAL") {
    datos_clima$year  <- as.numeric(substr(as.character(datos_clima$date), 1, 4))
    datos_clima$month <- as.numeric(substr(as.character(datos_clima$date), 5, 6))
    env_data <- dcast(datos_clima, year ~ month, value.var = "valor_analisis")
    for (m in 1:12) if (!as.character(m) %in% colnames(env_data)) env_data[[as.character(m)]] <- NA
    env_data <- env_data[, c("year", as.character(1:12))]
    
  } else {
    datos_clima$date_fmt <- as.Date(as.character(datos_clima$date), format = "%Y%m%d")
    datos_clima$year <- format(datos_clima$date_fmt, "%Y")
    datos_clima$day_of_year <- as.numeric(format(datos_clima$date_fmt, "%j"))
    env_data <- dcast(datos_clima, year ~ day_of_year, value.var = "valor_analisis")
    max_days <- 366
    if (ncol(env_data) < max_days + 1) {
      env_data <- cbind(env_data, matrix(NA, nrow = nrow(env_data), ncol = max_days - ncol(env_data) + 1))
    }
    colnames(env_data) <- c("year", paste0("X", 1:max_days))
  }
  
  # 2. BUCLE CORRELACIONES CON FILTRO
  nombre_clima_base <- tools::file_path_sans_ext(basename(INPUT_CLIMA_FILE))
  
  archivos_cronos <- list.files(dir_cronos_output, pattern = "\\.txt$", full.names = TRUE)
  archivos_cronos <- archivos_cronos[grep(nombre_rwl_base, archivos_cronos)]
  
  if(length(archivos_cronos) == 0) stop("No hay cronologías para analizar.")
  
  dir.create(file.path(dir_corr_output, "txt"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir_corr_output, "plots"), recursive = TRUE, showWarnings = FALSE)
  
  for (ruta_crono in archivos_cronos) {
    nombre_archivo_solo <- tools::file_path_sans_ext(basename(ruta_crono))
    metodo_detectado <- sub(paste0(nombre_rwl_base, "_"), "", nombre_archivo_solo)
    
    # --- FILTRO DE INTERRUPTORES ---
    estado_metodo <- USAR_METODOS[[metodo_detectado]]
    if (is.null(estado_metodo) || isFALSE(estado_metodo)) next 
    
    # Construir ID de Análisis
    id_analisis <- paste0(nombre_rwl_base, "_", metodo_detectado, "_vs_", nombre_clima_base)
    
    # --- AÑADIR AÑOS AL NOMBRE SI HAY FILTRO ---
    if (FILTRAR_ANIOS) {
      id_analisis <- paste0(id_analisis, "_", ANIO_INICIO, "-", ANIO_FIN)
    }
    
    message("    --> Analizando: ", id_analisis)
    
    response <- read.table(ruta_crono, header = TRUE, sep = "\t")
    if(!"year" %in% names(response) && "Year" %in% names(response)) names(response)[names(response)=="Year"] <- "year"
    
    if ("res" %in% names(response) && !all(is.na(response$res))) {
      response <- response %>% select(year, res); names(response)[2] <- "TRW"; tipo_serie <- "RES"
    } else {
      response <- response %>% select(year, std); names(response)[2] <- "TRW"; tipo_serie <- "STD"
    }
    
    # --- ALINEACIÓN CON FILTRO DE AÑOS ---
    common_years <- intersect(response$year, env_data$year)
    
    if (FILTRAR_ANIOS) {
      target_years <- ANIO_INICIO:ANIO_FIN
      common_years <- intersect(common_years, target_years)
      if (length(common_years) == 0) {
        message("        [SKIP] El rango de años seleccionado no tiene datos comunes.")
        next
      }
      message("        >>> Periodo filtrado: ", min(common_years), " - ", max(common_years))
    }
    
    if(length(common_years) < 15) { message("        [SKIP] Años insuficientes (<15)."); next }
    
    response_final <- response[response$year %in% common_years, ]
    env_data_final <- env_data[env_data$year %in% common_years, ]
    row.names(response_final) <- response_final$year; response_final <- response_final %>% select(-year)
    row.names(env_data_final) <- env_data_final$year; env_data_final <- env_data_final %>% select(-year)
    
    tryCatch({
      if (TIPO_DATOS == "DIARIO") {
        resultados <- daily_response(
          alpha = ALPHA_LEVEL, response = response_final, env_data = env_data_final,
          method = "cor", metric = "r.squared", cor_method = "pearson", 
          lower_limit = 1, upper_limit = 732,
          row_names_subset = TRUE, remove_insignificant = TRUE, 
          previous_year = TRUE, reference_window = "end", aggregate_function = "mean"
        )
      } else {
        resultados <- monthly_response(
          alpha = ALPHA_LEVEL, response = response_final, env_data = env_data_final,
          method = "cor", metric = "r.squared", cor_method = "pearson", 
          lower_limit = 1, upper_limit = 12,
          row_names_subset = TRUE, remove_insignificant = TRUE, 
          previous_year = TRUE, reference_window = "end", aggregate_function = "mean"
        )
      }
      
      sufijo <- paste0(id_analisis, "_", tipo_serie)
      
      write.table(as.data.frame(resultados$calculations), 
                  file.path(dir_corr_output, "txt", paste0("Correl_", sufijo, ".txt")), 
                  sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
      
      colors <- colourvalues::colour_values(1:100, palette = "diverge_hsv")
      resultados$plot_heatmap <- resultados$plot_heatmap + scale_fill_gradientn(colors = colors) + 
        ggtitle(paste0(id_analisis, " [", tipo_serie, "]"))
      
      pdf(file.path(dir_corr_output, "plots", paste0("Heatmap_", sufijo, ".pdf")), width = 10, height = 8)
      plot(resultados$plot_heatmap); dev.off()
      
      pdf(file.path(dir_corr_output, "plots", paste0("Extreme_", sufijo, ".pdf")), width = 10, height = 6)
      plot(resultados$plot_extreme); dev.off()
      
    }, error = function(e) { message("        [ERROR] Fallo: ", e$message) })
  }
  
  message(">>> PARTE 2 COMPLETADA.")
  message("    Resultados en: ", dir_corr_output)
}

message("\n=== SCRIPT MAESTRO FINALIZADO ===")
