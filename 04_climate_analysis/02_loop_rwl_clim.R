# ==============================================================================
# 17_loop_rwl_clim.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Sequential Multi-Method Standardisation & Correlation Loop.
#   A comprehensive script that iterates through every raw RWL file, applies multiple 
#   standardisation (detrending) methods, and correlates each resulting index 
#   series with climate data.
#
#   Workflow per RWL File:
#   1. Load Raw RWL.
#   2. Apply Standardisation Methods:
#      - Adaptive Spline (Age-dependent stiffness).
#      - Fixed Spline (32 years).
#      - ModNegExp (Modified Negative Exponential).
#      - Mean Detrending.
#   3. Generate Chronologies (Standard/Residual) for each method.
#   4. Correlate with Climate:
#      - Uses 'dendroTools::daily_response' for moving window analysis.
#      - Generates heatmaps of daily correlations.
#
#   Use Case:
#   - Best for detailed, step-by-step analysis or debugging. For large batches, 
#     use the parallel version (18_loop_rwl_clim_parallel.R).
#
#   Inputs:
#   - Directory of Raw RWL files.
#   - Directory of Climate Data files.
#
#   Outputs:
#   - Chronologies (.txt, .pdf, .rwi) for each method.
#   - Correlation Results (Heatmaps, Text files).
# ==============================================================================


# --- PARÁMETROS A MODIFICAR ---
# --- PARÁMETROS A MODIFICAR ---
INPUT_DIR_RWL       <- 'PLACEHOLDER/path/to/input_rwl_dir'
INPUT_DIR_CLIM      <- 'PLACEHOLDER/path/to/input_clim_dir'
OUTPUT_DIR_BASE     <- 'PLACEHOLDER/path/to/output_base_dir'

ALPHA_CORRELATION   <- 0.05

# Métodos de estandarización
STANDARDIZATION_METHODS <- list(
  list(id = "AdapSpl", name = "Adaptive Spline (original)", type = "adaptive"),
  list(id = "Spl32",   name = "Spline 32yr (fixed)",    type = "dplR", dplR_method = "Spline", dplR_nyrs = 32),
  list(id = "Spl67P",  name = "Spline 67% N",           type = "dplR", dplR_method = "Spline", dplR_nyrs_prop = 0.67),
  list(id = "ModNeg",  name = "ModNegExp",              type = "dplR", dplR_method = "ModNegExp"),
  list(id = "Mean",    name = "Mean Detrend",           type = "dplR", dplR_method = "Mean")
)
# ------------------------------


# --- (1) INSTALAR Y CARGAR LIBRERÍAS ---

# --- (1) INSTALAR Y CARGAR LIBRERÍAS ---

library(dplR)
library(dendroTools)
library(ggplot2)
library(dplyr)
library(reshape2)
library(colourvalues)
library(tools)

# Crear directorio de salida base si no existe
if (!dir.exists(OUTPUT_DIR_BASE)) {
  dir.create(OUTPUT_DIR_BASE, recursive = TRUE)
  message(paste("Directorio base de resultados creado:", OUTPUT_DIR_BASE))
}

# --- (PRE-LISTADO DE ARCHIVOS) ---
rwl_files_list <- list.files(path = INPUT_DIR_RWL, pattern = "\\.rwl$", full.names = TRUE)
if (length(rwl_files_list) == 0) {
  stop(paste("No se encontraron archivos .rwl en el directorio:", INPUT_DIR_RWL))
}
message(paste("\nSe encontraron", length(rwl_files_list), "archivos .rwl para procesar."))

clim_files_list <- list.files(path = INPUT_DIR_CLIM, pattern = "\\.txt$", full.names = TRUE)
if (length(clim_files_list) == 0) {
  stop(paste("No se encontraron archivos climáticos .txt en el directorio:", INPUT_DIR_CLIM))
}
message(paste("Se encontraron", length(clim_files_list), "archivos climáticos .txt para correlacionar."))


# --- (2) BUCLE PRINCIPAL PARA PROCESAR CADA ARCHIVO RWL ---
for (current_file_path_rwl in rwl_files_list) {
  
  message(paste("\n\n================ PROCESANDO ARCHIVO RWL:", basename(current_file_path_rwl), "================"))
  rwl_file_root <- sub("\\.rwl$", "", basename(current_file_path_rwl))
  
  # Leer datos RWL una vez por archivo
  raw_rwl_data <- NULL
  real_years_rwl <- NULL # Años reales del archivo RWL
  temp_df_for_detrending_original <- NULL # DataFrame original para detrending
  
  tryCatch({
    raw_rwl_data <- read.rwl(current_file_path_rwl)
    message(paste("Datos RWL cargados desde:", current_file_path_rwl))
    
    real_years_rwl <- as.numeric(rownames(raw_rwl_data))
    if(length(real_years_rwl) == 0 || all(is.na(real_years_rwl))) {
      stop("real_years_rwl no se pudo generar correctamente. Verifica los rownames del archivo RWL.")
    }
    # Preparar el data.frame para detrending (se usará en el bucle de estandarizaciones)
    temp_df_for_detrending_original <- as.data.frame(raw_rwl_data)
    rownames(temp_df_for_detrending_original) <- real_years_rwl
    
  }, error = function(e) {
    message(paste("\nERROR FATAL al cargar o pre-procesar el archivo RWL:", basename(current_file_path_rwl)))
    message(paste("Mensaje de error:", e$message))
    message("Saltando TODAS las estandarizaciones y correlaciones para este archivo RWL.\n")
    next # Salta al siguiente archivo RWL en el bucle exterior
  })
  
  if (is.null(raw_rwl_data)) next # Si falló la carga, ya se imprimió mensaje, saltar.
  
  # --- (2.1) BUCLE PARA CADA MÉTODO DE ESTANDARIZACIÓN ---
  for (std_method_params in STANDARDIZATION_METHODS) {
    
    std_id_short <- std_method_params$id
    std_name_descriptive <- std_method_params$name
    
    message(paste("\n  --- Aplicando estandarización:", std_name_descriptive, "(ID:", std_id_short, ") para", rwl_file_root, "---"))
    
    # Nombres de archivo y prefijos específicos para esta estandarización
    rwl_file_root_with_std_id <- paste0(rwl_file_root, "_", std_id_short)
    
    output_chronology_txt_current_std <- NULL # Para almacenar la ruta del archivo de cronología generado para ESTA estandarización
    chronology_generated_successfully_this_std <- FALSE
    
    tryCatch({
      # --- (2.1.1) DETRENDING Y GENERACIÓN DE CRONOLOGÍA PARA EL MÉTODO ACTUAL ---
      message(paste("  --- Iniciando Parte 1: Generación de Cronología (", std_id_short, ") para", rwl_file_root, "---"))
      
      output_dir_chron <- file.path(OUTPUT_DIR_BASE, "Chronologies")
      if (!dir.exists(output_dir_chron)) dir.create(output_dir_chron, recursive = TRUE)
      
      output_chronology_txt_current_std <- file.path(output_dir_chron, paste0(rwl_file_root_with_std_id, "_chronology.txt"))
      output_chronology_pdf_current_std <- file.path(output_dir_chron, paste0(rwl_file_root_with_std_id, "_chronology_plot.pdf"))
      output_rwi_filepath_current_std <- file.path(output_dir_chron, paste0(rwl_file_root_with_std_id, ".rwi"))
      
      # Función adaptativa original (definida fuera para que sea accesible)
      detrend_adaptativo_original <- function(series_vector, series_id, all_years_func_arg) {
        series_clean <- series_vector[!is.na(series_vector)]
        series_length <- length(series_clean)
        
        if (is.null(names(series_vector))) names(series_vector) <- all_years_func_arg
        
        if (series_length == 0) {
          message(paste("Advertencia: La serie", series_id, "está vacía. Devolviendo NAs."))
          empty_series <- rep(NA_real_, length(all_years_func_arg)); names(empty_series) <- all_years_func_arg
          return(empty_series)
        }
        
        if (series_length >= 100) detrended <- detrend.series(y = series_vector, method = "Spline", nyrs = 32, return.info = FALSE)
        else if (series_length >= 50) detrended <- detrend.series(y = series_vector, method = "Spline", nyrs = 20, return.info = FALSE)
        else if (series_length >= 30) {
          detrended <- tryCatch(detrend.series(y = series_vector, method = "ModNegExp", return.info = FALSE),
                                error = function(e) {
                                  message(paste("ModNegExp falló (", series_id, "): usando Spline nyrs=10. E:", e$message))
                                  detrend.series(y = series_vector, method = "Spline", nyrs = 10, return.info = FALSE)
                                })
        } else {
          detrended <- tryCatch(detrend.series(y = series_vector, method = "ModNegExp", return.info = FALSE),
                                error = function(e) {
                                  message(paste("ModNegExp falló (", series_id, "): usando Spline nyrs=5 o Mean. E:", e$message))
                                  if (series_length < 5) {
                                    message(paste("Serie", series_id, "<5, usando Mean."))
                                    detrend.series(y = series_vector, method = "Mean", return.info = FALSE)
                                  } else detrend.series(y = series_vector, method = "Spline", nyrs = 5, return.info = FALSE)
                                })
        }
        return(detrended)
      }
      
      detrended_series_list <- NULL
      if (std_method_params$type == "adaptive") {
        detrended_series_list <- lapply(seq_along(temp_df_for_detrending_original), function(i) {
          series_vector <- temp_df_for_detrending_original[[i]]
          series_id_map <- colnames(temp_df_for_detrending_original)[i]
          names(series_vector) <- rownames(temp_df_for_detrending_original) # Asegurar nombres de años
          detrend_adaptativo_original(series_vector, series_id = series_id_map, all_years_func_arg = rownames(temp_df_for_detrending_original))
        })
      } else if (std_method_params$type == "dplR") {
        detrended_series_list <- lapply(seq_along(temp_df_for_detrending_original), function(i) {
          series_vector <- temp_df_for_detrending_original[[i]]
          series_id_map <- colnames(temp_df_for_detrending_original)[i]
          names(series_vector) <- rownames(temp_df_for_detrending_original) # Asegurar nombres de años
          
          series_clean <- series_vector[!is.na(series_vector)]
          if (length(series_clean) == 0) {
            message(paste("Advertencia: Serie", series_id_map, "vacía para método", std_id_short,". Devolviendo NAs."))
            empty_series <- rep(NA_real_, length(real_years_rwl)); names(empty_series) <- real_years_rwl
            return(empty_series)
          }
          
          args_detrend <- list(y = series_vector, method = std_method_params$dplR_method, return.info = FALSE)
          
          if (!is.null(std_method_params$dplR_nyrs)) {
            args_detrend$nyrs <- std_method_params$dplR_nyrs
          } else if (!is.null(std_method_params$dplR_nyrs_prop)) {
            series_length <- length(series_clean)
            calculated_nyrs <- max(5, round(series_length * std_method_params$dplR_nyrs_prop)) # min 5 años para spline
            args_detrend$nyrs <- calculated_nyrs
            # message(paste("  Serie:", series_id_map, "Long:", series_length, "Calculated nyrs for SplP:", calculated_nyrs)) # Debug
          }
          
          # Para métodos como "Mean" o "ModNegExp", dplR::detrend.series no usa/necesita 'nyrs' explícitamente o lo ignora.
          # Solo pasamos 'nyrs' si el método es "Spline".
          if (std_method_params$dplR_method != "Spline" && "nyrs" %in% names(args_detrend)) {
            args_detrend$nyrs <- NULL # Eliminar nyrs si no es Spline
          }
          
          # Asegurar que para series muy cortas, ciertos métodos no fallen o se use uno más robusto
          if (length(series_clean) < 5 && std_method_params$dplR_method %in% c("Spline", "ModNegExp")) {
            message(paste("Serie", series_id_map, "demasiado corta (<5) para", std_method_params$dplR_method, ", usando detrend por media para esta serie."))
            args_detrend$method <- "Mean"
            if ("nyrs" %in% names(args_detrend)) args_detrend$nyrs <- NULL
          }
          
          ds_out <- tryCatch(
            do.call(dplR::detrend.series, args_detrend),
            error = function(e_ds) {
              message(paste("Error detrending serie", series_id_map, "con método", std_id_short, "(", args_detrend$method, "):", e_ds$message))
              message("  Intentando con detrend por Media como fallback para esta serie.")
              tryCatch(
                dplR::detrend.series(y = series_vector, method = "Mean", return.info = FALSE),
                error = function(e_mean) {
                  message(paste("Fallback a Mean también falló para serie", series_id_map, ":", e_mean$message, "- Devolviendo NAs."))
                  empty_series <- rep(NA_real_, length(real_years_rwl)); names(empty_series) <- real_years_rwl
                  return(empty_series)
                }
              )
            }
          )
          return(ds_out)
        })
      }
      names(detrended_series_list) <- colnames(temp_df_for_detrending_original)
      
      detrended_data_current_std <- NULL
      tryCatch({
        detrended_data_current_std <- as.data.frame(detrended_series_list, row.names = real_years_rwl)
      }, error = function(e){
        message(paste("  Error al convertir lista de series a data.frame:", e$message))
        message("  Intentando merge más robusto...")
        all_years_df_std <- data.frame(year = real_years_rwl)
        merged_data_std <- all_years_df_std
        for (s_name_std in names(detrended_series_list)) {
          if (!is.null(detrended_series_list[[s_name_std]])) { 
            series_df_std <- data.frame(year = as.numeric(names(detrended_series_list[[s_name_std]])),
                                        value = detrended_series_list[[s_name_std]])
            colnames(series_df_std)[2] <- s_name_std
            merged_data_std <- merge(merged_data_std, series_df_std, by = "year", all.x = TRUE)
          } else {
            message(paste("  Serie", s_name_std, "es NULL y será omitida del merge para", std_id_short))
          }
        }
        rownames(merged_data_std) <- merged_data_std$year
        detrended_data_current_std <<- merged_data_std[, -1, drop = FALSE]
      })
      
      if(is.null(detrended_data_current_std) || nrow(detrended_data_current_std) == 0 || ncol(detrended_data_current_std) == 0) {
        stop(paste("  Falló la creación de detrended_data para", rwl_file_root_with_std_id, "o resultó en un data.frame vacío."))
      }
      message(paste("  Detrending (", std_id_short, ") completado para todas las series."))
      
      chronology_object_current_std <- chron(detrended_data_current_std, prewhiten = TRUE, prefix = rwl_file_root_with_std_id) # Usar prefijo con std_id
      message(paste("  Cronologías (", std_id_short, ") estándar y residual generadas."))
      
      save_chronology_custom <- function(chronology_obj, file_path_chr, root_name_prefix_chr) {
        chronology_df_save <- as.data.frame(chronology_obj)
        # Los nombres de columna generados por chron() usan el prefijo dado
        std_col_name_prefixed_save <- paste0(root_name_prefix_chr, "std") 
        res_col_name_prefixed_save <- paste0(root_name_prefix_chr, "res")
        
        df_to_save_chr <- data.frame(year = as.numeric(rownames(chronology_df_save)))
        
        if (std_col_name_prefixed_save %in% colnames(chronology_df_save)) {
          df_to_save_chr$std <- chronology_df_save[[std_col_name_prefixed_save]]
        } else if ("std" %in% colnames(chronology_df_save)) { # Fallback si el prefijo no funcionó como se esperaba
          df_to_save_chr$std <- chronology_df_save[["std"]]
          message("  Advertencia: Se usó 'std' genérico para cronología estándar (prefijo no encontrado).")
        } else {
          message(paste("  Advertencia: No se encontró columna de cronología estándar para", root_name_prefix_chr))
          df_to_save_chr$std <- NA 
        }
        
        if (res_col_name_prefixed_save %in% colnames(chronology_df_save)) {
          df_to_save_chr$res <- chronology_df_save[[res_col_name_prefixed_save]]
        } else if ("res" %in% colnames(chronology_df_save)) { # Fallback
          df_to_save_chr$res <- chronology_df_save[["res"]]
          message("  Advertencia: Se usó 'res' genérico para cronología residual (prefijo no encontrado).")
        } else {
          message(paste("  Advertencia: No se encontró columna de cronología residual para", root_name_prefix_chr))
          df_to_save_chr$res <- NA 
        }
        
        if ("samp.depth" %in% colnames(chronology_df_save)) {
          df_to_save_chr$samp_depth <- chronology_df_save[["samp.depth"]]
        } else {
          message(paste("  Advertencia: No se encontró 'samp.depth' para", root_name_prefix_chr))
          df_to_save_chr$samp_depth <- NA 
        }
        
        write.table(df_to_save_chr, file = file_path_chr, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        message(paste("  Cronologías (", std_id_short, ") guardadas en:", file_path_chr))
      }
      
      save_chronology_custom(chronology_object_current_std, output_chronology_txt_current_std, rwl_file_root_with_std_id)
      
      pdf(output_chronology_pdf_current_std, width = 14, height = 8)
      plot(chronology_object_current_std, add.spline = TRUE, nyrs = 20,
           main = paste("Cronología (", std_name_descriptive, ") para", rwl_file_root),
           xlab = "Años", ylab = "Índice de Ancho de Anillo (RWI)")
      dev.off()
      message(paste("  Gráfico de cronología (", std_id_short, ") guardado en:", output_chronology_pdf_current_std))
      
      colnames(detrended_data_current_std) <- gsub("\\.dat$", "", colnames(detrended_data_current_std)) # Limpiar nombres si es necesario
      write.rwl(as.data.frame(detrended_data_current_std), fname = output_rwi_filepath_current_std, format = "tucson")
      message(paste("  Archivo RWI (", std_id_short, ") guardado en:", output_rwi_filepath_current_std))
      
      chronology_generated_successfully_this_std <- TRUE
      message(paste("  --- Fin Parte 1: Generación de Cronología (", std_id_short, ") ---"))
      
    }, error = function(e) {
      message(paste("\n  ERROR en Parte 1 (Generación de Cronología con método", std_id_short, ") para:", rwl_file_root))
      message(paste("  Mensaje de error:", e$message))
      message(paste("  No se podrá realizar la correlación para esta combinación RWL-Estandarización (",std_id_short,").\n"))
      chronology_generated_successfully_this_std <- FALSE
    }) # Fin de tryCatch para Parte 1 (generación de crono para esta estandarización)
    
    # --- (2.1.2) BUCLE PARA PROCESAR CADA ARCHIVO CLIMÁTICO CONTRA LA CRONOLOGÍA RWL-STD ACTUAL ---
    if (chronology_generated_successfully_this_std && !is.null(output_chronology_txt_current_std)) {
      
      for (current_file_path_clim in clim_files_list) {
        
        clim_file_basename <- basename(current_file_path_clim)
        clim_variable_name_short <- sub("([^_.]+).*", "\\1", clim_file_basename)
        
        message(paste("\n    --- Iniciando Correlación de", rwl_file_root_with_std_id, "con archivo climático:", clim_file_basename, "---"))
        
        tryCatch({
          # PRE-PROCESAMIENTO DE DATOS CLIMÁTICOS (se repite por si acaso, pero podría optimizarse si siempre es el mismo)
          # Esto es rápido, así que no es un gran cuello de botella repetirlo.
          clim_data_daily_raw <- read.table(current_file_path_clim, header = TRUE, sep = "\t",
                                            stringsAsFactors = FALSE, na.strings = c("NA", "N/A", "-", "", "null", "NaN"))
          if (ncol(clim_data_daily_raw) < 2) stop(paste("Archivo clima",clim_file_basename,"no tiene >= 2 columnas."))
          clim_value_col_name <- colnames(clim_data_daily_raw)[2]
          if (!is.numeric(clim_data_daily_raw[[clim_value_col_name]])) {
            clim_data_daily_raw[[clim_value_col_name]] <- suppressWarnings(as.numeric(as.character(clim_data_daily_raw[[clim_value_col_name]])))
            if(all(is.na(clim_data_daily_raw[[clim_value_col_name]]))) stop(paste("Falló conversión a numérico de '", clim_value_col_name, "' en", clim_file_basename))
          }
          clim_data_daily_filtered <- clim_data_daily_raw[!is.na(clim_data_daily_raw[[clim_value_col_name]]), ]
          colnames(clim_data_daily_filtered)[1] <- "date_str"
          clim_data_daily_filtered$date <- as.Date(as.character(clim_data_daily_filtered$date_str), format = "%Y%m%d")
          if(any(is.na(clim_data_daily_filtered$date))) stop(paste("Fechas en '", clim_file_basename, "' no parseadas (YYYYMMDD)."))
          clim_data_daily_filtered <- clim_data_daily_filtered[!is.na(clim_data_daily_filtered$date), ]
          clim_data_daily_filtered$year <- format(clim_data_daily_filtered$date, "%Y")
          clim_data_daily_filtered$day_of_year_X <- paste0("X", as.numeric(format(clim_data_daily_filtered$date, "%j")))
          env_data_clim_prepared <- dcast(clim_data_daily_filtered, year ~ day_of_year_X, value.var = clim_value_col_name)
          max_days <- 366; all_doy_cols_X <- paste0("X", 1:max_days)
          for(col_name_X in all_doy_cols_X) if(!(col_name_X %in% colnames(env_data_clim_prepared))) env_data_clim_prepared[[col_name_X]] <- NA_real_
          env_data_clim_prepared <- env_data_clim_prepared[, c("year", all_doy_cols_X)]
          message(paste("    Datos de", clim_variable_name_short, "cargados y formateados."))
          
          # PARTE 2: CORRELACIÓN
          chron_data_for_corr_std <- read.table(output_chronology_txt_current_std, header = TRUE, sep = "\t")
          if (!"year" %in% names(chron_data_for_corr_std)) stop("Columna 'year' no encontrada en cronología.")
          names(chron_data_for_corr_std)[names(chron_data_for_corr_std) == "year"] <- "Year"
          
          if (!"res" %in% names(chron_data_for_corr_std) || all(is.na(chron_data_for_corr_std$res))) {
            message(paste("    ADVERTENCIA: Col 'res' en", basename(output_chronology_txt_current_std), "es NA. Saltando correlación."))
            next 
          }
          response_chron_std <- chron_data_for_corr_std %>% select(Year, res) %>% filter(!is.na(res))
          if (nrow(response_chron_std) == 0) {
            message(paste("    ADVERTENCIA: No hay datos 'res' válidos en", basename(output_chronology_txt_current_std), ". Saltando correlación."))
            next 
          }
          if (!is.numeric(response_chron_std$res)) stop("Col 'res' no es numérica.")
          
          names(response_chron_std)[names(response_chron_std) == "res"] <- "TRW"
          response_chron_std$Year <- as.numeric(response_chron_std$Year)
          
          common_years_std <- intersect(as.character(response_chron_std$Year), env_data_clim_prepared$year)
          if(length(common_years_std) < 3) stop(paste("Menos de 3 años comunes entre crono", rwl_file_root_with_std_id, "y clima", clim_variable_name_short))
          
          response_chron_aligned_std <- response_chron_std[as.character(response_chron_std$Year) %in% common_years_std, ]
          env_data_clim_aligned_current_std <- env_data_clim_prepared[env_data_clim_prepared$year %in% common_years_std, ] 
          
          rownames(response_chron_aligned_std) <- response_chron_aligned_std$Year
          response_chron_aligned_std <- response_chron_aligned_std %>% select(TRW)
          rownames(env_data_clim_aligned_current_std) <- env_data_clim_aligned_current_std$year
          env_data_clim_aligned_current_std <- env_data_clim_aligned_current_std %>% select(-year)
          env_data_clim_aligned_current_std <- as.data.frame(sapply(env_data_clim_aligned_current_std, as.numeric))
          rownames(env_data_clim_aligned_current_std) <- common_years_std
          message("    Datos de cronología y clima alineados.")
          
          resultados_corr_std <- daily_response(
            response = response_chron_aligned_std,
            env_data = env_data_clim_aligned_current_std,
            method = "cor", metric = "r.squared", lower_limit = 1, upper_limit = 732,
            cor_method = "pearson", alpha = ALPHA_CORRELATION, row_names_subset = TRUE,
            remove_insignificant = TRUE, previous_year = TRUE, reference_window = "end", aggregate_function = "mean"
          )
          message("    Análisis de correlación diario completado.")
          
          # Rutas de salida para correlaciones, ahora con std_id
          ruta_salida_corr_base_clim_std <- file.path(OUTPUT_DIR_BASE, paste0("Correlacion_", clim_variable_name_short)) # Carpeta por variable climática
          ruta_salida_corr_txt_clim_std <- file.path(ruta_salida_corr_base_clim_std, "txt")
          ruta_salida_corr_plots_clim_std <- file.path(ruta_salida_corr_base_clim_std, "plots")
          if (!dir.exists(ruta_salida_corr_txt_clim_std)) dir.create(ruta_salida_corr_txt_clim_std, recursive = TRUE)
          if (!dir.exists(ruta_salida_corr_plots_clim_std)) dir.create(ruta_salida_corr_plots_clim_std, recursive = TRUE)
          
          alpha_str_corr <- sprintf("alfa%02d", ALPHA_CORRELATION * 100)
          output_filename_base_corr_std <- paste0(rwl_file_root_with_std_id, "_", clim_variable_name_short, "_correl_", alpha_str_corr) # Incluye std_id
          
          archivo_resultados_txt_std <- file.path(ruta_salida_corr_txt_clim_std, paste0(output_filename_base_corr_std, ".txt"))
          if (!is.null(resultados_corr_std$calculations) && nrow(resultados_corr_std$calculations) > 0) {
            write.table(as.data.frame(resultados_corr_std$calculations), archivo_resultados_txt_std, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
            message(paste("    Resultados de correlación guardados en:", archivo_resultados_txt_std))
          } else {
            cat(paste("No se encontraron resultados significativos para guardar."), file = archivo_resultados_txt_std)
            message(paste("    No se encontraron resultados significativos para:", archivo_resultados_txt_std))
          }
          
          if (!is.null(resultados_corr_std$plot_heatmap)) {
            custom_heatmap_std <- resultados_corr_std$plot_heatmap +
              scale_fill_gradientn(colours = colourvalues::colour_values(1:100, palette = "diverge_hsv"), na.value = "grey90") +
              ggtitle(paste("Correlación Diaria:", rwl_file_root, "(", std_id_short, ") vs", clim_variable_name_short)) # Título con std_id
            
            ggsave(filename = file.path(ruta_salida_corr_plots_clim_std, paste0(output_filename_base_corr_std, "_heatmap.pdf")), plot = custom_heatmap_std, width = 12, height = 7)
            ggsave(filename = file.path(ruta_salida_corr_plots_clim_std, paste0(output_filename_base_corr_std, "_heatmap.png")), plot = custom_heatmap_std, width = 12, height = 7, dpi = 300)
            message(paste("    Gráfico heatmap (",std_id_short,") guardado."))
          }
          
          if (!is.null(resultados_corr_std$plot_extreme) && inherits(resultados_corr_std$plot_extreme, "ggplot")) {
            # Modificar título del plot_extreme si es necesario (dendroTools puede no permitirlo directamente)
            # Podrías intentar: resultados_corr_std$plot_extreme <- resultados_corr_std$plot_extreme + ggtitle(...) pero podría no funcionar con todos los plots.
            ggsave(filename = file.path(ruta_salida_corr_plots_clim_std, paste0(output_filename_base_corr_std, "_extreme.pdf")), plot = resultados_corr_std$plot_extreme, width = 10, height = 6)
            ggsave(filename = file.path(ruta_salida_corr_plots_clim_std, paste0(output_filename_base_corr_std, "_extreme.png")), plot = resultados_corr_std$plot_extreme, width = 10, height = 6, dpi=300)
            message(paste("    Gráfico de extremos (",std_id_short,") guardado."))
          }
          
          message(paste0("    --- Fin Correlación de ", rwl_file_root_with_std_id, " con ", clim_variable_name_short, " ---"))
          
        }, error = function(e) {
          message(paste("\n    ERROR procesando correlación para climático:", clim_file_basename, "con RWL-Std:", rwl_file_root_with_std_id))
          message(paste("    Mensaje de error:", e$message))
          message("    Saltando a la siguiente combinación.\n")
        }) # Fin de tryCatch para el archivo climático actual
      } # Fin del bucle FOR sobre archivos climáticos
    } else {
      message(paste("  Saltando correlaciones para", rwl_file_root_with_std_id, "debido a error en su generación de cronología."))
    } # Fin de if (chronology_generated_successfully_this_std)
  } # Fin del bucle FOR sobre métodos de estandarización
  
  message(paste("================ FIN PROCESAMIENTO RWL:", basename(current_file_path_rwl), "================"))
} # Fin del bucle FOR sobre archivos rwl

message("\n--- Script completado para todos los archivos .rwl, métodos de estandarización y variables climáticas ---")
