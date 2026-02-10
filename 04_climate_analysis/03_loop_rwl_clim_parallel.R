# ==============================================================================
# 18_loop_rwl_clim_parallel.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Parallel Multi-Method Standardisation & Correlation Loop.
#   The high-performance version of '17_loop_rwl_clim.R', designed to process large 
#   datasets efficiently by distributing work across multiple CPU cores.
#
#   Features:
#   - Uses the `future` and `future.apply` packages for parallel execution.
#   - Replicates the full workflow of the sequential script:
#     (RWL Loading -> Multi-Method Detrending -> Chronology Building -> Climate Correlation).
#   - Robust error handling within each worker to prevent full crash on single file failure.
#
#   Configuration:
#   - Automatically detects available cores (leaves 1 free for OS).
#   - NUM_WORKERS parameter allows manual control.
#
#   Inputs/Outputs:
#   - Identical to 17_loop_rwl_clim.R
# ==============================================================================


# Cargar librerías
library(dplR)
library(dendroTools)
library(ggplot2)
library(dplyr)
library(reshape2)
library(colourvalues)
library(tools)
library(future)
library(future.apply)
library(parallel)

# --- PARÁMETROS A MODIFICAR ---
INPUT_DIR_RWL       <- 'PLACEHOLDER/path/to/input_rwl_dir'
INPUT_DIR_CLIM      <- 'PLACEHOLDER/path/to/input_clim_dir'
OUTPUT_DIR_BASE     <- 'PLACEHOLDER/path/to/output_base_dir'

ALPHA_CORRELATION   <- 0.05
NUM_WORKERS         <- max(1, parallel::detectCores() - 1)

# Métodos de estandarización
STANDARDIZATION_METHODS <- list(
  list(id = "AdapSpl", name = "Adaptive Spline (original)", type = "adaptive"),
  list(id = "Spl32",   name = "Spline 32yr (fixed)",    type = "dplR", dplR_method = "Spline", dplR_nyrs = 32),
  list(id = "Spl67P",  name = "Spline 67% N",           type = "dplR", dplR_method = "Spline", dplR_nyrs_prop = 0.67),
  list(id = "ModNeg",  name = "ModNegExp",              type = "dplR", dplR_method = "ModNegExp"),
  list(id = "Mean",    name = "Mean Detrend",           type = "dplR", dplR_method = "Mean")
process_single_rwl_file <- function(current_file_path_rwl, 
                                    standardisation_methods_list_arg, 
                                    clim_files_list_arg, 
                                    output_dir_base_arg, 
                                    alpha_correlation_arg,
                                    loaded_libraries_arg) {
  
  # Cargar librerías necesarias DENTRO de la función para los workers paralelos
  for (lib_worker in loaded_libraries_arg) {
    suppressPackageStartupMessages(library(lib_worker, character.only = TRUE, quietly = TRUE))
  }
  
  worker_pid <- Sys.getpid()
  message(paste0("[Worker PID: ", worker_pid, "] PROCESANDO ARCHIVO RWL: ", basename(current_file_path_rwl)))
  rwl_file_root <- sub("\\.rwl$", "", basename(current_file_path_rwl))
  
  raw_rwl_data <- NULL
  real_years_rwl <- NULL
  temp_df_for_detrending_original <- NULL
  
  # Definición de la función adaptativa DENTRO del worker para asegurar disponibilidad
  detrend_adaptativo_original <- function(series_vector, series_id, all_years_func_arg) {
    series_clean <- series_vector[!is.na(series_vector)]
    series_length <- length(series_clean)
    
    if (is.null(names(series_vector))) names(series_vector) <- all_years_func_arg
    
    if (series_length == 0) {
      message(paste0("  [Worker PID:", worker_pid, "] Advertencia: La serie ", series_id, " está vacía. Devolviendo NAs."))
      empty_series <- rep(NA_real_, length(all_years_func_arg)); names(empty_series) <- all_years_func_arg
      return(empty_series)
    }
    
    if (series_length >= 100) detrended <- detrend.series(y = series_vector, method = "Spline", nyrs = 32, return.info = FALSE)
    else if (series_length >= 50) detrended <- detrend.series(y = series_vector, method = "Spline", nyrs = 20, return.info = FALSE)
    else if (series_length >= 30) {
      detrended <- tryCatch(detrend.series(y = series_vector, method = "ModNegExp", return.info = FALSE),
                            error = function(e) {
                              message(paste0("  [Worker PID:", worker_pid, "] ModNegExp falló (", series_id, "): usando Spline nyrs=10. E:", e$message))
                              detrend.series(y = series_vector, method = "Spline", nyrs = 10, return.info = FALSE)
                            })
    } else {
      detrended <- tryCatch(detrend.series(y = series_vector, method = "ModNegExp", return.info = FALSE),
                            error = function(e) {
                              message(paste0("  [Worker PID:", worker_pid, "] ModNegExp falló (", series_id, "): usando Spline nyrs=5 o Mean. E:", e$message))
                              if (series_length < 5) {
                                message(paste0("  [Worker PID:", worker_pid, "] Serie ", series_id, "<5, usando Mean."))
                                detrend.series(y = series_vector, method = "Mean", return.info = FALSE)
                              } else detrend.series(y = series_vector, method = "Spline", nyrs = 5, return.info = FALSE)
                            })
    }
    return(detrended)
  }
  
  # Función para guardar cronología, definida una vez dentro del worker
  save_chronology_custom_worker <- function(chronology_obj, file_path_chr, root_name_prefix_chr, std_id_chr) {
    chronology_df_save <- as.data.frame(chronology_obj)
    std_col_name_prefixed_save <- paste0(root_name_prefix_chr, "std") 
    res_col_name_prefixed_save <- paste0(root_name_prefix_chr, "res")
    
    df_to_save_chr <- data.frame(year = as.numeric(rownames(chronology_df_save)))
    
    if (std_col_name_prefixed_save %in% colnames(chronology_df_save)) {
      df_to_save_chr$std <- chronology_df_save[[std_col_name_prefixed_save]]
    } else if ("std" %in% colnames(chronology_df_save)) {
      df_to_save_chr$std <- chronology_df_save[["std"]]
      message(paste0("  [Worker PID:", worker_pid, "] Advertencia (",std_id_chr,"): Se usó 'std' genérico para cronología estándar (prefijo no encontrado)."))
    } else {
      message(paste0("  [Worker PID:", worker_pid, "] Advertencia (",std_id_chr,"): No se encontró columna de cronología estándar para ", root_name_prefix_chr))
      df_to_save_chr$std <- NA 
    }
    
    if (res_col_name_prefixed_save %in% colnames(chronology_df_save)) {
      df_to_save_chr$res <- chronology_df_save[[res_col_name_prefixed_save]]
    } else if ("res" %in% colnames(chronology_df_save)) {
      df_to_save_chr$res <- chronology_df_save[["res"]]
      message(paste0("  [Worker PID:", worker_pid, "] Advertencia (",std_id_chr,"): Se usó 'res' genérico para cronología residual (prefijo no encontrado)."))
    } else {
      message(paste0("  [Worker PID:", worker_pid, "] Advertencia (",std_id_chr,"): No se encontró columna de cronología residual para ", root_name_prefix_chr))
      df_to_save_chr$res <- NA 
    }
    
    if ("samp.depth" %in% colnames(chronology_df_save)) {
      df_to_save_chr$samp_depth <- chronology_df_save[["samp.depth"]]
    } else {
      message(paste0("  [Worker PID:", worker_pid, "] Advertencia (",std_id_chr,"): No se encontró 'samp.depth' para ", root_name_prefix_chr))
      df_to_save_chr$samp_depth <- NA 
    }
    
    write.table(df_to_save_chr, file = file_path_chr, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    message(paste0("  [Worker PID:", worker_pid, "] Cronologías (", std_id_chr, ") guardadas en: ", file_path_chr))
  }
  
  
  tryCatch({
    raw_rwl_data <- read.rwl(current_file_path_rwl)
    message(paste0("  [Worker PID:", worker_pid, "] Datos RWL cargados desde: ", current_file_path_rwl))
    
    real_years_rwl <- as.numeric(rownames(raw_rwl_data))
    if(length(real_years_rwl) == 0 || all(is.na(real_years_rwl))) {
      stop("real_years_rwl no se pudo generar correctamente. Verifica los rownames del archivo RWL.")
    }
    temp_df_for_detrending_original <- as.data.frame(raw_rwl_data)
    rownames(temp_df_for_detrending_original) <- real_years_rwl
    
  }, error = function(e) {
    message(paste0("\n  [Worker PID:", worker_pid, "] ERROR FATAL al cargar o pre-procesar el archivo RWL: ", basename(current_file_path_rwl)))
    message(paste0("  [Worker PID:", worker_pid, "] Mensaje de error: ", e$message))
    message(paste0("  [Worker PID:", worker_pid, "] Saltando TODAS las estandarizaciones y correlaciones para este archivo RWL.\n"))
    return(list(file = basename(current_file_path_rwl), status = "ERROR_RWL_LOAD", message = e$message))
  })
  
  if (is.null(raw_rwl_data)) {
    return(list(file = basename(current_file_path_rwl), status = "SKIPPED_RWL_LOAD_FAILED", message = "raw_rwl_data is NULL after tryCatch."))
  }
  
  # --- (2.1) BUCLE PARA CADA MÉTODO DE ESTANDARIZACIÓN ---
  for (std_method_params in STANDARDIZATION_METHODS) {
    
    std_id_short <- std_method_params$id
    std_name_descriptive <- std_method_params$name
    
    message(paste0("\n  [Worker PID:", worker_pid, "] --- Aplicando estandarización: ", std_name_descriptive, " (ID: ", std_id_short, ") para ", rwl_file_root, " ---"))
    
    rwl_file_root_with_std_id <- paste0(rwl_file_root, "_", std_id_short)
    
    output_chronology_txt_current_std <- NULL
    chronology_generated_successfully_this_std <- FALSE
    
    # Variable para el dataframe de series destendenciadas (definida aquí para que esté en el scope del <<- si se usa)
    detrended_data_current_std <- NULL 
    
    tryCatch({
      message(paste0("  [Worker PID:", worker_pid, "] --- Iniciando Parte 1: Generación de Cronología (", std_id_short, ") para ", rwl_file_root, " ---"))
      
      output_dir_chron <- file.path(OUTPUT_DIR_BASE, "Chronologies")
      if (!dir.exists(output_dir_chron)) dir.create(output_dir_chron, recursive = TRUE)
      
      output_chronology_txt_current_std <- file.path(output_dir_chron, paste0(rwl_file_root_with_std_id, "_chronology.txt"))
      output_chronology_pdf_current_std <- file.path(output_dir_chron, paste0(rwl_file_root_with_std_id, "_chronology_plot.pdf"))
      output_rwi_filepath_current_std <- file.path(output_dir_chron, paste0(rwl_file_root_with_std_id, ".rwi"))
      
      detrended_series_list <- NULL
      if (std_method_params$type == "adaptive") {
        detrended_series_list <- lapply(seq_along(temp_df_for_detrending_original), function(i) {
          series_vector <- temp_df_for_detrending_original[[i]]
          series_id_map <- colnames(temp_df_for_detrending_original)[i]
          names(series_vector) <- rownames(temp_df_for_detrending_original)
          detrend_adaptativo_original(series_vector, series_id = series_id_map, all_years_func_arg = rownames(temp_df_for_detrending_original))
        })
      } else if (std_method_params$type == "dplR") {
        detrended_series_list <- lapply(seq_along(temp_df_for_detrending_original), function(i) {
          series_vector <- temp_df_for_detrending_original[[i]]
          series_id_map <- colnames(temp_df_for_detrending_original)[i]
          names(series_vector) <- rownames(temp_df_for_detrending_original)
          
          series_clean <- series_vector[!is.na(series_vector)]
          if (length(series_clean) == 0) {
            message(paste0("  [Worker PID:", worker_pid, "] Advertencia: Serie ", series_id_map, " vacía para método ", std_id_short,". Devolviendo NAs."))
            empty_series <- rep(NA_real_, length(real_years_rwl)); names(empty_series) <- real_years_rwl
            return(empty_series)
          }
          
          args_detrend <- list(y = series_vector, method = std_method_params$dplR_method, return.info = FALSE)
          
          if (!is.null(std_method_params$dplR_nyrs)) {
            args_detrend$nyrs <- std_method_params$dplR_nyrs
          } else if (!is.null(std_method_params$dplR_nyrs_prop)) {
            series_length <- length(series_clean)
            calculated_nyrs <- max(5, round(series_length * std_method_params$dplR_nyrs_prop))
            args_detrend$nyrs <- calculated_nyrs
          }
          
          if (std_method_params$dplR_method != "Spline" && "nyrs" %in% names(args_detrend)) {
            args_detrend$nyrs <- NULL
          }
          
          if (length(series_clean) < 5 && std_method_params$dplR_method %in% c("Spline", "ModNegExp")) {
            message(paste0("  [Worker PID:", worker_pid, "] Serie ", series_id_map, " demasiado corta (<5) para ", std_method_params$dplR_method, ", usando detrend por media para esta serie."))
            args_detrend$method <- "Mean"
            if ("nyrs" %in% names(args_detrend)) args_detrend$nyrs <- NULL
          }
          
          ds_out <- tryCatch(
            do.call(dplR::detrend.series, args_detrend),
            error = function(e_ds) {
              message(paste0("  [Worker PID:", worker_pid, "] Error detrending serie ", series_id_map, " con método ", std_id_short, " (", args_detrend$method, "): ", e_ds$message))
              message(paste0("  [Worker PID:", worker_pid, "]   Intentando con detrend por Media como fallback para esta serie."))
              tryCatch(
                dplR::detrend.series(y = series_vector, method = "Mean", return.info = FALSE),
                error = function(e_mean) {
                  message(paste0("  [Worker PID:", worker_pid, "]   Fallback a Mean también falló para serie ", series_id_map, ": ", e_mean$message, " - Devolviendo NAs."))
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
      
      # detrended_data_current_std ya está definida en el scope de esta función.
      tryCatch({
        detrended_data_current_std <- as.data.frame(detrended_series_list, row.names = real_years_rwl)
      }, error = function(e){
        message(paste0("  [Worker PID:", worker_pid, "] Error al convertir lista de series a data.frame (",std_id_short,"): ", e$message))
        message(paste0("  [Worker PID:", worker_pid, "]   Intentando merge más robusto (",std_id_short,")..."))
        all_years_df_std <- data.frame(year = real_years_rwl)
        merged_data_std <- all_years_df_std
        for (s_name_std in names(detrended_series_list)) {
          if (!is.null(detrended_series_list[[s_name_std]])) { 
            series_df_std <- data.frame(year = as.numeric(names(detrended_series_list[[s_name_std]])),
                                        value = detrended_series_list[[s_name_std]])
            colnames(series_df_std)[2] <- s_name_std
            merged_data_std <- merge(merged_data_std, series_df_std, by = "year", all.x = TRUE)
          } else {
            message(paste0("  [Worker PID:", worker_pid, "]   Serie ", s_name_std, " es NULL y será omitida del merge para ", std_id_short))
          }
        }
        rownames(merged_data_std) <- merged_data_std$year
        detrended_data_current_std <<- merged_data_std[, -1, drop = FALSE] # Asigna a la variable en el scope de este bucle de estandarización
      })
      
      if(is.null(detrended_data_current_std) || nrow(detrended_data_current_std) == 0 || ncol(detrended_data_current_std) == 0) {
        stop(paste0("Falló la creación de detrended_data para ", rwl_file_root_with_std_id, " o resultó en un data.frame vacío."))
      }
      message(paste0("  [Worker PID:", worker_pid, "] Detrending (", std_id_short, ") completado para todas las series."))
      
      chronology_object_current_std <- chron(detrended_data_current_std, prewhiten = TRUE, prefix = rwl_file_root_with_std_id)
      message(paste0("  [Worker PID:", worker_pid, "] Cronologías (", std_id_short, ") estándar y residual generadas."))
      
      save_chronology_custom_worker(chronology_object_current_std, output_chronology_txt_current_std, rwl_file_root_with_std_id, std_id_short)
      
      pdf(output_chronology_pdf_current_std, width = 14, height = 8)
      plot(chronology_object_current_std, add.spline = TRUE, nyrs = 20,
           main = paste("Cronología (", std_name_descriptive, ") para", rwl_file_root),
           xlab = "Años", ylab = "Índice de Ancho de Anillo (RWI)")
      dev.off()
      message(paste0("  [Worker PID:", worker_pid, "] Gráfico de cronología (", std_id_short, ") guardado en: ", output_chronology_pdf_current_std))
      
      colnames(detrended_data_current_std) <- gsub("\\.dat$", "", colnames(detrended_data_current_std))
      write.rwl(as.data.frame(detrended_data_current_std), fname = output_rwi_filepath_current_std, format = "tucson")
      message(paste0("  [Worker PID:", worker_pid, "] Archivo RWI (", std_id_short, ") guardado en: ", output_rwi_filepath_current_std))
      
      chronology_generated_successfully_this_std <- TRUE
      message(paste0("  [Worker PID:", worker_pid, "] --- Fin Parte 1: Generación de Cronología (", std_id_short, ") ---"))
      
    }, error = function(e) {
      message(paste0("\n  [Worker PID:", worker_pid, "] ERROR en Parte 1 (Generación de Cronología con método ", std_id_short, ") para: ", rwl_file_root))
      message(paste0("  [Worker PID:", worker_pid, "] Mensaje de error: ", e$message))
      message(paste0("  [Worker PID:", worker_pid, "] No se podrá realizar la correlación para esta combinación RWL-Estandarización (",std_id_short,").\n"))
      chronology_generated_successfully_this_std <- FALSE
    })
    
    # --- (2.1.2) BUCLE PARA PROCESAR CADA ARCHIVO CLIMÁTICO CONTRA LA CRONOLOGÍA RWL-STD ACTUAL ---
    if (chronology_generated_successfully_this_std && !is.null(output_chronology_txt_current_std)) {
      
      for (current_file_path_clim in clim_files_list_arg) {
        
        clim_file_basename <- basename(current_file_path_clim)
        clim_variable_name_short <- sub("([^_.]+).*", "\\1", clim_file_basename)
        
        message(paste0("\n    [Worker PID:", worker_pid, "] --- Iniciando Correlación de ", rwl_file_root_with_std_id, " con archivo climático: ", clim_file_basename, " ---"))
        
        tryCatch({
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
          message(paste0("    [Worker PID:", worker_pid, "] Datos de ", clim_variable_name_short, " cargados y formateados."))
          
          chron_data_for_corr_std <- read.table(output_chronology_txt_current_std, header = TRUE, sep = "\t")
          if (!"year" %in% names(chron_data_for_corr_std)) stop("Columna 'year' no encontrada en cronología.")
          names(chron_data_for_corr_std)[names(chron_data_for_corr_std) == "year"] <- "Year"
          
          if (!"res" %in% names(chron_data_for_corr_std) || all(is.na(chron_data_for_corr_std$res))) {
            message(paste0("    [Worker PID:", worker_pid, "] ADVERTENCIA: Col 'res' en ", basename(output_chronology_txt_current_std), " es NA. Saltando correlación."))
            next 
          }
          response_chron_std <- chron_data_for_corr_std %>% select(Year, res) %>% filter(!is.na(res))
          if (nrow(response_chron_std) == 0) {
            message(paste0("    [Worker PID:", worker_pid, "] ADVERTENCIA: No hay datos 'res' válidos en ", basename(output_chronology_txt_current_std), ". Saltando correlación."))
            next 
          }
          if (!is.numeric(response_chron_std$res)) stop("Col 'res' no es numérica.")
          
          names(response_chron_std)[names(response_chron_std) == "res"] <- "TRW"
          response_chron_std$Year <- as.numeric(response_chron_std$Year)
          
          common_years_std <- intersect(as.character(response_chron_std$Year), env_data_clim_prepared$year)
          if(length(common_years_std) < 3) stop(paste0("Menos de 3 años comunes entre crono ", rwl_file_root_with_std_id, " y clima ", clim_variable_name_short))
          
          response_chron_aligned_std <- response_chron_std[as.character(response_chron_std$Year) %in% common_years_std, ]
          env_data_clim_aligned_current_std <- env_data_clim_prepared[env_data_clim_prepared$year %in% common_years_std, ] 
          
          rownames(response_chron_aligned_std) <- response_chron_aligned_std$Year
          response_chron_aligned_std <- response_chron_aligned_std %>% select(TRW)
          rownames(env_data_clim_aligned_current_std) <- env_data_clim_aligned_current_std$year
          env_data_clim_aligned_current_std <- env_data_clim_aligned_current_std %>% select(-year)
          env_data_clim_aligned_current_std <- as.data.frame(sapply(env_data_clim_aligned_current_std, as.numeric))
          rownames(env_data_clim_aligned_current_std) <- common_years_std
          message(paste0("    [Worker PID:", worker_pid, "] Datos de cronología y clima alineados."))
          
          resultados_corr_std <- daily_response(
            response = response_chron_aligned_std,
            env_data = env_data_clim_aligned_current_std,
            method = "cor", metric = "r.squared", lower_limit = 1, upper_limit = 732,
            cor_method = "pearson", alpha = ALPHA_CORRELATION, row_names_subset = TRUE,
            remove_insignificant = TRUE, previous_year = TRUE, reference_window = "end", aggregate_function = "mean"
          )
          message(paste0("    [Worker PID:", worker_pid, "] Análisis de correlación diario completado."))
          
          ruta_salida_corr_base_clim_std <- file.path(OUTPUT_DIR_BASE, paste0("Correlacion_", clim_variable_name_short))
          ruta_salida_corr_txt_clim_std <- file.path(ruta_salida_corr_base_clim_std, "txt")
          ruta_salida_corr_plots_clim_std <- file.path(ruta_salida_corr_base_clim_std, "plots")
          if (!dir.exists(ruta_salida_corr_txt_clim_std)) dir.create(ruta_salida_corr_txt_clim_std, recursive = TRUE)
          if (!dir.exists(ruta_salida_corr_plots_clim_std)) dir.create(ruta_salida_corr_plots_clim_std, recursive = TRUE)
          
          alpha_str_corr <- sprintf("alfa%02d", ALPHA_CORRELATION * 100)
          output_filename_base_corr_std <- paste0(rwl_file_root_with_std_id, "_", clim_variable_name_short, "_correl_", alpha_str_corr)
          
          archivo_resultados_txt_std <- file.path(ruta_salida_corr_txt_clim_std, paste0(output_filename_base_corr_std, ".txt"))
          if (!is.null(resultados_corr_std$calculations) && nrow(resultados_corr_std$calculations) > 0) {
            write.table(as.data.frame(resultados_corr_std$calculations), archivo_resultados_txt_std, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
            message(paste0("    [Worker PID:", worker_pid, "] Resultados de correlación guardados en: ", archivo_resultados_txt_std))
          } else {
            cat(paste("No se encontraron resultados significativos para guardar."), file = archivo_resultados_txt_std)
            message(paste0("    [Worker PID:", worker_pid, "] No se encontraron resultados significativos para: ", archivo_resultados_txt_std))
          }
          
          if (!is.null(resultados_corr_std$plot_heatmap)) {
            custom_heatmap_std <- resultados_corr_std$plot_heatmap +
              scale_fill_gradientn(colours = colourvalues::colour_values(1:100, palette = "diverge_hsv"), na.value = "grey90") +
              ggtitle(paste("Correlación Diaria:", rwl_file_root, "(", std_id_short, ") vs", clim_variable_name_short))
            
            ggsave(filename = file.path(ruta_salida_corr_plots_clim_std, paste0(output_filename_base_corr_std, "_heatmap.pdf")), plot = custom_heatmap_std, width = 12, height = 7)
            ggsave(filename = file.path(ruta_salida_corr_plots_clim_std, paste0(output_filename_base_corr_std, "_heatmap.png")), plot = custom_heatmap_std, width = 12, height = 7, dpi = 300)
            message(paste0("    [Worker PID:", worker_pid, "] Gráfico heatmap (",std_id_short,") guardado."))
          }
          
          if (!is.null(resultados_corr_std$plot_extreme) && inherits(resultados_corr_std$plot_extreme, "ggplot")) {
            ggsave(filename = file.path(ruta_salida_corr_plots_clim_std, paste0(output_filename_base_corr_std, "_extreme.pdf")), plot = resultados_corr_std$plot_extreme, width = 10, height = 6)
            ggsave(filename = file.path(ruta_salida_corr_plots_clim_std, paste0(output_filename_base_corr_std, "_extreme.png")), plot = resultados_corr_std$plot_extreme, width = 10, height = 6, dpi=300)
            message(paste0("    [Worker PID:", worker_pid, "] Gráfico de extremos (",std_id_short,") guardado."))
          }
          
          message(paste0("    [Worker PID:", worker_pid, "] --- Fin Correlación de ", rwl_file_root_with_std_id, " con ", clim_variable_name_short, " ---"))
          
        }, error = function(e) {
          message(paste0("\n    [Worker PID:", worker_pid, "] ERROR procesando correlación para climático: ", clim_file_basename, " con RWL-Std: ", rwl_file_root_with_std_id))
          message(paste0("    [Worker PID:", worker_pid, "] Mensaje de error: ", e$message))
          message(paste0("    [Worker PID:", worker_pid, "] Saltando a la siguiente combinación.\n"))
        })
      }
    } else {
      message(paste0("  [Worker PID:", worker_pid, "] Saltando correlaciones para ", rwl_file_root_with_std_id, " debido a error en su generación de cronología."))
    }
  }
  
  message(paste0("[Worker PID:", worker_pid, "] ================ FIN PROCESAMIENTO RWL: ", basename(current_file_path_rwl), " ================"))
  return(list(file = basename(current_file_path_rwl), status = "SUCCESS", message = "Procesado exitosamente."))
}


# --- (2) EJECUTAR PROCESAMIENTO EN PARALELO ---
# Configurar el plan de `future`.
# `multisession` crea workers R en segundo plano.
# `plan(sequential)` puede usarse para depurar sin paralelismo.
if (NUM_WORKERS > 1) {
  plan(multisession, workers = NUM_WORKERS)
  message(paste("\nIniciando procesamiento paralelo con", future::nbrOfWorkers(), "workers."))
} else {
  plan(sequential)
  message("\nIniciando procesamiento en modo secuencial (1 worker).")
}


# Usar future_lapply para aplicar la función a cada archivo RWL
# future.seed = TRUE es importante para la reproducibilidad.
# Es buena práctica pasar explícitamente todas las variables que la función necesita.
processing_results <- future_lapply(
  X = rwl_files_list, 
  FUN = process_single_rwl_file,
  standardisation_methods_list_arg = STANDARDIZATION_METHODS,
  clim_files_list_arg = clim_files_list,
  output_dir_base_arg = OUTPUT_DIR_BASE,
  alpha_correlation_arg = ALPHA_CORRELATION,
  loaded_libraries_arg = libraries_to_load, # Pasar la lista de librerías
  future.seed = TRUE # Para reproducibilidad
  # No es necesario future.globals si las funciones como detrend_adaptativo_original
  # están definidas DENTRO de process_single_rwl_file o son parte de los paquetes cargados.
)

message("\n--- Procesamiento paralelo (o secuencial con estructura future) completado ---")

# (Opcional) Inspeccionar los resultados/errores devueltos por cada worker
message("\nResumen de resultados del procesamiento:")
successful_files <- 0
failed_files <- 0
for (i in seq_along(processing_results)) {
  result_item <- processing_results[[i]]
  message(paste("Archivo:", result_item$file, "- Estado:", result_item$status))
  if (grepl("ERROR|FAIL|SKIP", result_item$status)) {
    failed_files <- failed_files + 1
    message(paste("  Detalle:", result_item$message))
  } else {
    successful_files <- successful_files + 1
  }
}
message(paste("\nTotal archivos procesados con éxito:", successful_files))
message(paste("Total archivos con fallos o saltados:", failed_files))

message("\n--- Script completado para todos los archivos .rwl, métodos de estandarización y variables climáticas ---")

# Es buena práctica volver al modo secuencial, lo que cierra los workers de fondo.
plan(sequential)
