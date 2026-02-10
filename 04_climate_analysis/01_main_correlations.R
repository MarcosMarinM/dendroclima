# ==============================================================================
# 16_main_correlations.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Dendroclimatic Correlation Analysis (Core Script).
#   The primary analysis engine for calculating relationships between tree-ring 
#   chronologies and climate variables.
#
#   Functionality:
#   1. Data Loading:
#      - Reads formatted daily climate data (converted to monthly/seasonal).
#      - Loads all processed tree-ring chronologies (.txt format).
#   2. Correlation Analysis:
#      - Computes Pearson correlations for:
#        - Monthly windows (Previous Year Oct-Dec, Current Year Jan-Dec).
#        - Seasonal aggregates (Winter, Spring, Summer, Autumn).
#        - Annual averages.
#   3. Visualisation:
#      - Generates detailed static line plots for every climate variable (Precip, Tmax, Tmin, etc.)
#        and every chronology, split by Age Group (Young, Mid, Old).
#   4. SPEI Integration:
#      - Specifically handles SPEI (Standardised Precipitation Evapotranspiration Index) 
#        correlations across multiple time scales (1, 3, 6, 12, etc.).
#
#   Inputs:
#   - Directory of Climate Text Files.
#   - Directory of Chronology Text Files.
#   - SPEI Directory.
#
#   Outputs:
#   - Correlation Matrices (.rds).
#   - PDF Plots of Correlations.
# ==============================================================================

# --- PARÁMETROS A MODIFICAR ---
# 0. CARGA DE LIBRERÍAS
library(dplyr)
library(tidyr)
library(ggplot2)
library(colourvalues)
library(stringr)

# --- PARÁMETROS A MODIFICAR ---
PATH_CLIM_DIR          <- "PLACEHOLDER/path/to/climate_dir/" 
PATH_CHRON_TXT_DIR     <- "PLACEHOLDER/path/to/chronologies/"
OUTPUT_DIR_CLIMATE_CORR <- 'PLACEHOLDER/path/to/output_correlations'
PATH_SPEI_DIR          <- "PLACEHOLDER/path/to/spei_dir"

# Crear directorio de salida si no existe
dir.create(OUTPUT_DIR_CLIMATE_CORR, showWarnings = FALSE, recursive = TRUE)
# ------------------------------



# --- 1. Definición de Archivos Climáticos ---
files_clim <- list(
  precip = file.path(PATH_CLIM_DIR, "precip_data.txt"),
  tmax   = file.path(PATH_CLIM_DIR, "tmax_data.txt"),
  tmed   = file.path(PATH_CLIM_DIR, "tmed_data.txt"),
  tmin   = file.path(PATH_CLIM_DIR, "tmin_data.txt")
)

# Definition of chronology groups for explorations
chrono_sets <- list(
  set1 = c("TRW", "LW", "EW"),
  set2 = c("LWBI", "EWBI", "DBI")
)

# --- 2. Helper Functions ---
source("../00_Utils/00_utils_dendro.R")


# --- 3. Process Climate Data ---
cat("\n--- Processing Climate Data (original) ---\n")
climate_data_daily <- list()
# Note: The original script used "precipitation" for precip, "tmax" for tmax, etc.
# these are the names that will be stored in `climate_data_daily` and `climate_data_aggregated`.
climate_data_daily$precip <- read_daily_climate_data(files_clim$precip, "precipitation")
climate_data_daily$tmax   <- read_daily_climate_data(files_clim$tmax, "tmax")
climate_data_daily$tmed   <- read_daily_climate_data(files_clim$tmed, "tmed")
climate_data_daily$tmin   <- read_daily_climate_data(files_clim$tmin, "tmin")

climate_data_aggregated <- list()
climate_data_aggregated$precip <- aggregate_climate_data(climate_data_daily$precip, agg_fun = "sum")
climate_data_aggregated$tmax   <- aggregate_climate_data(climate_data_daily$tmax, agg_fun = "mean")
climate_data_aggregated$tmed   <- aggregate_climate_data(climate_data_daily$tmed, agg_fun = "mean")
climate_data_aggregated$tmin   <- aggregate_climate_data(climate_data_daily$tmin, agg_fun = "mean")

# --- 4. Load and Prepare Tree-Ring Chronologies ---
cat("\n--- Loading Tree-Ring Chronologies ---\n")
chron_files_txt <- list.files(path = PATH_CHRON_TXT_DIR, pattern = "\\.txt$", full.names = TRUE)
if (length(chron_files_txt) == 0) stop("No .txt chronologies found.")

all_chronologies <- list()
for (chron_file in chron_files_txt) {
  chron_name_only <- basename(chron_file)
  info <- get_info_from_chron_filename(chron_name_only)
  if (info$age_group %in% c("unknown","all") || info$type == "UNKNOWN") {
    warning(paste("Chronology", chron_name_only, "skipped (non-standard age/type)."))
    next
  }
  current_chron_data <- tryCatch(read.table(chron_file, header=TRUE, sep="\t", stringsAsFactors=FALSE),
                                 error = function(e) {warning(paste("Error reading", chron_name_only)); NULL})
  if (is.null(current_chron_data) || !all(c("year","res") %in% colnames(current_chron_data))) {
    warning(paste("Incorrect format in", chron_name_only, "- skipping.")); next
  }
  if (is.null(all_chronologies[[info$age_group]])) all_chronologies[[info$age_group]] <- list()
  all_chronologies[[info$age_group]][[info$type]] <- select(current_chron_data, year, res)
}
cat(length(unlist(all_chronologies, recursive=FALSE)), "chronologies loaded.\n")

# --- 5. Calculate Correlations (for original climate variables) ---
cat("\n--- Calculating Correlations (for original climate variables) ---\n")

# Define labels for correlation periods
prev_months_labels_full <- paste0("p", month.abb) # pJan, pFeb, ... (for internal calculation)
curr_months_labels_full <- month.abb             # Jan, Feb, ... (for internal calculation)
season_labels <- c("DJF", "MAM", "JJA", "SON")
annual_label <- "Annual"

correlation_results <- list()

for (age_g in names(all_chronologies)) {
  correlation_results[[age_g]] <- list()
  for (chrono_t in names(all_chronologies[[age_g]])) {
    correlation_results[[age_g]][[chrono_t]] <- list()
    current_chrono_series <- all_chronologies[[age_g]][[chrono_t]]
    if(is.null(current_chrono_series) || nrow(current_chrono_series) < 10) next
    
    for (clim_var_name_key in names(climate_data_aggregated)) {
      correlation_results[[age_g]][[chrono_t]][[clim_var_name_key]] <- list()
      clim_monthly_data <- climate_data_aggregated[[clim_var_name_key]]$monthly
      
      # Current year months
      for (m_idx in 1:12) {
        period_label <- curr_months_labels_full[m_idx]
        clim_s <- clim_monthly_data %>% filter(month == m_idx) %>% select(year, value)
        common <- inner_join(current_chrono_series, clim_s, by="year")
        r_val <- if(nrow(common) >= 10 && sd(common$res, na.rm=T)>0 && sd(common$value, na.rm=T)>0) cor(common$res, common$value, use="p", method="p") else NA
        correlation_results[[age_g]][[chrono_t]][[clim_var_name_key]][[period_label]] <- r_val
      }
      # Previous year months (all 12, but we'll plot only Oct-Dec)
      for (m_idx in 1:12) {
        period_label <- prev_months_labels_full[m_idx] 
        clim_s <- clim_monthly_data %>% filter(month == m_idx) %>% mutate(year = year + 1) %>% select(year, value)
        common <- inner_join(current_chrono_series, clim_s, by="year")
        r_val <- if(nrow(common) >= 10 && sd(common$res, na.rm=T)>0 && sd(common$value, na.rm=T)>0) cor(common$res, common$value, use="p", method="p") else NA
        correlation_results[[age_g]][[chrono_t]][[clim_var_name_key]][[period_label]] <- r_val
      }
      
      # Seasonal (current year)
      clim_seasonal_data <- climate_data_aggregated[[clim_var_name_key]]$seasonal 
      for (s_lab in season_labels) {
        clim_s <- clim_seasonal_data %>% filter(season == s_lab) %>% select(year, value)
        common <- inner_join(current_chrono_series, clim_s, by="year")
        r_val <- if(nrow(common) >= 10 && sd(common$res, na.rm=T)>0 && sd(common$value, na.rm=T)>0) cor(common$res, common$value, use="p", method="p") else NA
        correlation_results[[age_g]][[chrono_t]][[clim_var_name_key]][[s_lab]] <- r_val
      }
      
      # Annual (current year)
      clim_annual_data <- climate_data_aggregated[[clim_var_name_key]]$annual 
      common <- inner_join(current_chrono_series, clim_annual_data, by="year")
      r_val <- if(nrow(common) >= 10 && sd(common$res, na.rm=T)>0 && sd(common$value, na.rm=T)>0) cor(common$res, common$value, use="p", method="p") else NA
      correlation_results[[age_g]][[chrono_t]][[clim_var_name_key]][[annual_label]] <- r_val
      
      cat("  Corr:", age_g, "-", chrono_t, "vs", clim_var_name_key, "OK\n")
    }
  }
}
cat("\n--- Correlation calculation completed. --- \n")

results_rds_file <- file.path(OUTPUT_DIR_CLIMATE_CORR, "correlation_results_full.rds")
saveRDS(correlation_results, file = results_rds_file)
cat("Results saved to:", results_rds_file, "\n")


# --- 6. Generate Specific Correlation Plots (for original climate variables) ---
cat("\n--- Generating Correlation Plots by Age Group (original climate vars) ---\n")

# Load results if not in memory (useful if running script in chunks)
if (!exists("correlation_results")) {
  if (file.exists(results_rds_file)) {
    cat("Loading results from:", results_rds_file, "\n")
    correlation_results <- readRDS(results_rds_file)
  } else {
    stop(paste("File", results_rds_file, "not found. Run the first part of the script."))
  }
}

# Prepare full data.frame for plotting (all climate variables)
plot_data_list <- list()
for (age_g in names(correlation_results)) {
  for (chrono_t in names(correlation_results[[age_g]])) {
    for (clim_var_key in names(correlation_results[[age_g]][[chrono_t]])) {
      clim_var_data <- correlation_results[[age_g]][[chrono_t]][[clim_var_key]]
      if (length(clim_var_data) > 0) {
        temp_df <- data.frame(
          age_group = age_g,
          chrono_type = chrono_t,
          climate_variable = clim_var_key,
          period = names(clim_var_data),
          correlation = unlist(clim_var_data),
          stringsAsFactors = FALSE
        )
        plot_data_list[[length(plot_data_list) + 1]] <- temp_df
      }
    }
  }
}

if (length(plot_data_list) == 0) {
  stop("No correlation data found to plot. Check 'correlation_results'.")
}
plot_data_df_all_clim <- bind_rows(plot_data_list)
rownames(plot_data_df_all_clim) <- NULL

valid_chronos <- unique(unlist(chrono_sets)) # Ensure only desired chronos are plotted
plot_data_df_all_clim <- plot_data_df_all_clim %>%
  filter(chrono_type %in% valid_chronos)

# Define order for periods to ensure correct x-axis order
# Only display pOct, pNov, pDec for previous year months
prev_months_labels_plot <- paste0("p", month.abb[10:12]) 
curr_months_labels_plot <- month.abb             
season_labels_plot <- c("DJF", "MAM", "JJA", "SON")
annual_label_plot <- "Annual"
all_period_levels <- c(prev_months_labels_plot, curr_months_labels_plot, season_labels_plot, annual_label_plot)

plot_data_df_all_clim$period <- factor(plot_data_df_all_clim$period, levels = all_period_levels, ordered = TRUE)
plot_data_df_all_clim <- plot_data_df_all_clim %>% filter(!is.na(period)) # Remove any periods not in the defined levels

if (nrow(plot_data_df_all_clim) == 0) {
  stop("No data remains after filtering periods. Check period labels and data content.")
}

# Define exploration set for faceting
plot_data_df_all_clim <- plot_data_df_all_clim %>%
  mutate(
    exploration_set = case_when(
      chrono_type %in% chrono_sets$set1 ~ "Set 1: TRW, LW, EW",
      chrono_type %in% chrono_sets$set2 ~ "Set 2: LWBI, EWBI, DBI",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(exploration_set)) # Remove rows where exploration_set is NA

# Define order for chronology types and age groups for plotting consistency
chrono_type_levels_ordered <- c("TRW", "EW", "LW", "DBI", "EWBI", "LWBI")
plot_data_df_all_clim <- plot_data_df_all_clim %>%
  mutate(chrono_type = factor(chrono_type, levels = chrono_type_levels_ordered)) %>%
  filter(!is.na(chrono_type)) # Filter out any chrono_type not in the ordered list

age_group_levels <- c("yng", "mid", "old")
age_group_labels_en <- c("Young Trees", "Mid-Age Trees", "Old Trees")
if(all(unique(plot_data_df_all_clim$age_group) %in% age_group_levels)){
  plot_data_df_all_clim$age_group <- factor(plot_data_df_all_clim$age_group, levels = age_group_levels, labels = age_group_labels_en)
} else {
  cat("Warning: Age groups do not match 'yng', 'mid', 'old'. Original names will be used.\n")
  plot_data_df_all_clim$age_group <- factor(plot_data_df_all_clim$age_group)
}

# Custom colors and shapes for plotting
custom_colors <- c(
  "TRW" = "royalblue4", "EW" = "skyblue2", "LW" = "cadetblue1",
  "DBI" = "lightgreen", "EWBI" = "limegreen", "LWBI" = "darkgreen"
)
custom_colors <- custom_colors[chrono_type_levels_ordered] # Ensure order for legend

custom_shapes <- c(
  "TRW" = 16,  # Filled circle
  "EW"  = 15,  # Filled square
  "LW"  = 17,  # Filled triangle
  "DBI" = 16,  # Filled circle
  "EWBI"= 15,  # Filled square
  "LWBI"= 17   # Filled triangle
)
custom_shapes <- custom_shapes[chrono_type_levels_ordered] # Ensure order for legend

# Vertical line position to separate "previous year months" from "current year months"
vline_position_x_clim <- length(prev_months_labels_plot) + 0.5 # Line after pDec

background_color_aggregated <- "grey90" # For the aggregated periods (seasonal/annual)

climate_vars_to_plot <- c("precip", "tmax", "tmed", "tmin")
climate_vars_titles <- c(
  precip = "Precipitation",
  tmax   = "Maximum Temperature",
  tmed   = "Mean Temperature",
  tmin   = "Minimum Temperature"
)

for (current_clim_var_key in climate_vars_to_plot) {
  current_clim_var_title <- climate_vars_titles[current_clim_var_key]
  cat("\n--- Generating Plot for:", current_clim_var_title, "---\n")
  
  plot_data_df_current <- plot_data_df_all_clim %>%
    filter(climate_variable == current_clim_var_key)
  
  if (nrow(plot_data_df_current) == 0) {
    cat("No data found for", current_clim_var_title, "- skipping plot.\n")
    next
  }
  
  # Calculate dynamic Y-axis limits AND breaks for the current plot
  y_info <- get_dynamic_y_limits_and_breaks(plot_data_df_current$correlation)
  
  correlation_plot <- ggplot(plot_data_df_current, aes(x = period, y = correlation, group = chrono_type, color = chrono_type, shape = chrono_type)) +
    # Background for aggregated periods (seasonal & annual)
    geom_rect(
      xmin = length(prev_months_labels_plot) + length(curr_months_labels_plot) + 0.5, 
      xmax = Inf, ymin = -Inf, ymax = Inf,
      fill = background_color_aggregated, alpha = 0.3, color = NA, inherit.aes = FALSE
    ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(
      xintercept = vline_position_x_clim, linetype = "dotted", color = "black", linewidth = 0.8
    ) +
    facet_grid(age_group ~ exploration_set, scales = "free_y") +
    # Use dynamic limits and calculate breaks based on these limits
    scale_y_continuous(limits = y_info$limits, breaks = y_info$breaks) + 
    scale_color_manual(
      values = custom_colors,
      breaks = chrono_type_levels_ordered,
      name = "Wood Variable"
    ) +
    scale_shape_manual(
      values = custom_shapes,
      breaks = chrono_type_levels_ordered,
      name = "Wood Variable"
    ) +
    labs(
      title = paste("Correlation: Wood Variables and", current_clim_var_title),
      subtitle = "Monthly (prev. year Oct-Dec, curr. year Jan-Dec), Seasonal & Annual aggregates.",
      x = "Period",
      y = "Correlation Coefficient (r)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=7),
      axis.text.y = element_text(size=8),
      axis.title = element_text(size=10),
      plot.title = element_text(size=12, hjust = 0.5, face="bold"),
      plot.subtitle = element_text(size=9, hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size=9, face="bold"),
      legend.text = element_text(size=8),
      strip.text.x = element_text(size = 9, face="bold"),
      strip.text.y = element_text(size = 9, face="bold"),
      panel.spacing = unit(0.8, "lines")
    ) +
    guides(color = guide_legend(title="Wood Variable"), shape = guide_legend(title="Wood Variable"))
  
  pdf_filename_current <- paste0("Fig_Correlations_", current_clim_var_key, "_By_Age.pdf")
  output_pdf_path_current <- file.path(output_dir_climate_corr, pdf_filename_current)
  
  ggsave(
    filename = output_pdf_path_current,
    plot = correlation_plot,
    width = 11,
    height = 9,
    units = "in",
    device = "pdf"
  )
  cat("Plot for", current_clim_var_title, "saved to:", output_pdf_path_current, "\n")
}


# --- 7. NUEVA SECCIÓN: Procesamiento y Correlación de SPEI ---
cat("\n--- INICIANDO SECCIÓN DE SPEI ---\n")

# Parámetros específicos para SPEI
path_spei_dir <- "PLACEHOLDER/path/to/spei_dir" # !! AJUSTA ESTA RUTA !!
output_dir_spei_corr <- file.path(output_dir_climate_corr, "SPEI_Correlations") # Subcarpeta para SPEI
dir.create(output_dir_spei_corr, showWarnings = FALSE, recursive = TRUE)

# Define los meses del año previo y actual para los gráficos SPEI
# Estos son los meses que se considerarán para el eje X en el gráfico
prev_months_labels_spei_plot <- paste0("p", month.abb[10:12]) 
curr_months_labels_spei_plot <- month.abb             
all_period_levels_spei <- c(prev_months_labels_spei_plot, curr_months_labels_spei_plot)

# Función para leer archivos SPEI (sin encabezado de meses, solo año y luego meses)
read_spei_file <- function(file_path) {
  cat("  Reading SPEI file:", basename(file_path), "\n")
  # Nombres de columnas esperados: Año y los 12 meses
  col_names <- c("year", month.abb)
  
  spei_data <- tryCatch(read.table(file_path, header = FALSE, sep = "\t", 
                                   col.names = col_names, stringsAsFactors = FALSE, na.strings = c("NA", "-9999", "999.9", "-999.9")),
                        error = function(e) {
                          warning(paste("Error reading SPEI file", basename(file_path), ":", e$message))
                          return(NULL)
                        })
  
  if (is.null(spei_data) || nrow(spei_data) == 0) {
    warning(paste("No valid data found in SPEI file:", basename(file_path)))
    return(NULL)
  }
  
  # Pivotear de formato ancho a largo para facilitar joins
  spei_long <- spei_data %>%
    pivot_longer(
      cols = all_of(month.abb), # Pivotear las columnas de los meses
      names_to = "month_abb",
      values_to = "value"
    ) %>%
    mutate(month = match(month_abb, month.abb)) %>% # Convertir abreviatura de mes a número (1-12)
    select(year, month, value) %>%
    arrange(year, month)
  
  return(spei_long)
}

# Obtener lista de archivos SPEI
spei_files_full_path <- list.files(path = path_spei_dir, pattern = "^SPEI_\\d{2}meses\\.txt$", full.names = TRUE)

if (length(spei_files_full_path) == 0) {
  cat("ADVERTENCIA: No se encontraron archivos SPEI en la carpeta:", path_spei_dir, "\n")
  cat("  Saltando la sección de análisis y gráficos de SPEI.\n")
} else {
  cat("\n  Se encontraron", length(spei_files_full_path), "archivos SPEI para procesar.\n")
  
  all_spei_data_by_scale <- list()
  for (f_path in spei_files_full_path) {
    # Extraer la escala del nombre del archivo (e.g., "01" de "SPEI_01meses.txt")
    scale_match <- str_extract(basename(f_path), "(?<=SPEI_)\\d{2}(?=meses\\.txt)")
    if (!is.na(scale_match)) {
      scale_label <- paste0("SPEI_", as.integer(scale_match), "m")
      all_spei_data_by_scale[[scale_label]] <- read_spei_file(f_path)
    }
  }
  
  # Filtrar entradas NULL si hubo errores de lectura
  all_spei_data_by_scale <- Filter(Negate(is.null), all_spei_data_by_scale)
  
  if (length(all_spei_data_by_scale) == 0) {
    cat("ADVERTENCIA: Ningún archivo SPEI pudo ser leído correctamente. Saltando el análisis de SPEI.\n")
  } else {
    # Inicializar la lista de resultados de correlación para SPEI
    spei_correlation_results <- list()
    
    cat("\n--- Calculando correlaciones de SPEI ---\n")
    for (age_g in names(all_chronologies)) {
      spei_correlation_results[[age_g]] <- list()
      for (chrono_t in names(all_chronologies[[age_g]])) {
        spei_correlation_results[[age_g]][[chrono_t]] <- list()
        current_chrono_series <- all_chronologies[[age_g]][[chrono_t]]
        if(is.null(current_chrono_series) || nrow(current_chrono_series) < 10) next
        
        for (spei_scale_label in names(all_spei_data_by_scale)) {
          spei_correlation_results[[age_g]][[chrono_t]][[spei_scale_label]] <- list()
          current_spei_data_monthly <- all_spei_data_by_scale[[spei_scale_label]]
          
          # Correlacionar meses del AÑO ACTUAL
          for (m_idx in 1:12) {
            period_label <- month.abb[m_idx] # Jan, Feb, ...
            spei_s <- current_spei_data_monthly %>% filter(month == m_idx) %>% select(year, value)
            
            common <- inner_join(current_chrono_series, spei_s, by="year")
            r_val <- if(nrow(common) >= 10 && sd(common$res, na.rm=T)>0 && sd(common$value, na.rm=T)>0) cor(common$res, common$value, use="p", method="p") else NA
            
            spei_correlation_results[[age_g]][[chrono_t]][[spei_scale_label]][[period_label]] <- r_val
          }
          
          # Correlacionar meses del AÑO PREVIO (pOct, pNov, pDec)
          # Los meses 10, 11, 12 son Oct, Nov, Dec
          for (m_idx in 10:12) { 
            period_label <- paste0("p", month.abb[m_idx]) # pOct, pNov, pDec
            spei_s <- current_spei_data_monthly %>% 
              filter(month == m_idx) %>% 
              mutate(year = year + 1) %>% # Shift year for previous year's climate
              select(year, value)
            
            common <- inner_join(current_chrono_series, spei_s, by="year")
            r_val <- if(nrow(common) >= 10 && sd(common$res, na.rm=T)>0 && sd(common$value, na.rm=T)>0) cor(common$res, common$value, use="p", method="p") else NA
            
            spei_correlation_results[[age_g]][[chrono_t]][[spei_scale_label]][[period_label]] <- r_val
          }
          cat("    Corr:", age_g, "-", chrono_t, "vs", spei_scale_label, "OK\n")
        }
      }
    }
    
    # Guardar resultados de SPEI
    results_spei_rds_file <- file.path(output_dir_spei_corr, "correlation_results_SPEI.rds")
    saveRDS(spei_correlation_results, file = results_spei_rds_file)
    cat("Resultados de correlación SPEI guardados en:", results_spei_rds_file, "\n")
    
    # Preparar datos para los gráficos de SPEI
    plot_spei_data_list <- list()
    for (age_g in names(spei_correlation_results)) {
      for (chrono_t in names(spei_correlation_results[[age_g]])) {
        for (spei_scale_key in names(spei_correlation_results[[age_g]][[chrono_t]])) {
          spei_scale_data <- spei_correlation_results[[age_g]][[chrono_t]][[spei_scale_key]]
          if (length(spei_scale_data) > 0) {
            temp_df <- data.frame(
              age_group = age_g,
              chrono_type = chrono_t,
              spei_scale = spei_scale_key,
              period = names(spei_scale_data),
              correlation = unlist(spei_scale_data),
              stringsAsFactors = FALSE
            )
            plot_spei_data_list[[length(plot_spei_data_list) + 1]] <- temp_df
          }
        }
      }
    }
    
    if (length(plot_spei_data_list) == 0) {
      cat("ADVERTENCIA: No se encontraron datos de correlación SPEI para graficar. Saltando gráficos de SPEI.\n")
    } else {
      plot_data_df_spei <- bind_rows(plot_spei_data_list)
      rownames(plot_data_df_spei) <- NULL
      
      # Filtrar por cronologías válidas y definir sets de exploración
      plot_data_df_spei <- plot_data_df_spei %>%
        filter(chrono_type %in% valid_chronos) %>%
        mutate(
          exploration_set = case_when(
            chrono_type %in% chrono_sets$set1 ~ "Set 1: TRW, LW, EW",
            chrono_type %in% chrono_sets$set2 ~ "Set 2: LWBI, EWBI, DBI",
            TRUE ~ NA_character_
          ),
          chrono_type = factor(chrono_type, levels = chrono_type_levels_ordered)
        ) %>%
        filter(!is.na(exploration_set) & !is.na(chrono_type))
      
      # Aplicar niveles de grupo de edad
      if(all(unique(plot_data_df_spei$age_group) %in% age_group_levels)){
        plot_data_df_spei$age_group <- factor(plot_data_df_spei$age_group, levels = age_group_levels, labels = age_group_labels_en)
      } else {
        cat("Warning: Age groups for SPEI plots do not match 'yng', 'mid', 'old'. Original names will be used.\n")
        plot_data_df_spei$age_group <- factor(plot_data_df_spei$age_group)
      }
      
      # Definir el orden de los meses para el eje X
      plot_data_df_spei$period <- factor(plot_data_df_spei$period, levels = all_period_levels_spei, ordered = TRUE)
      plot_data_df_spei <- plot_data_df_spei %>% filter(!is.na(period))
      
      # Vertical line position to separate "previous year months" from "current year months" for SPEI plots
      vline_position_x_spei <- length(prev_months_labels_spei_plot) + 0.5 
      
      cat("\n--- Generando gráficos de correlación para SPEI ---\n")
      
      # Iterar sobre cada escala de SPEI para generar un gráfico separado
      for (current_spei_scale_label in unique(plot_data_df_spei$spei_scale)) {
        cat("  Generando gráfico para escala:", current_spei_scale_label, "\n")
        
        plot_data_df_current_spei_scale <- plot_data_df_spei %>%
          filter(spei_scale == current_spei_scale_label)
        
        if (nrow(plot_data_df_current_spei_scale) == 0) {
          cat("    No data found for SPEI scale", current_spei_scale_label, "- skipping plot.\n")
          next
        }
        
        # Calculate dynamic Y-axis limits AND breaks for the current SPEI plot
        y_info_spei <- get_dynamic_y_limits_and_breaks(plot_data_df_current_spei_scale$correlation)
        
        spei_correlation_plot <- ggplot(plot_data_df_current_spei_scale, aes(x = period, y = correlation, group = chrono_type, color = chrono_type, shape = chrono_type)) +
          geom_line(linewidth = 0.7) +
          geom_point(size = 1.6) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
          geom_vline(
            xintercept = vline_position_x_spei, linetype = "dotted", color = "black", linewidth = 0.8
          ) +
          facet_grid(age_group ~ exploration_set, scales = "free_y") +
          # Use dynamic limits and calculate breaks based on these limits
          scale_y_continuous(limits = y_info_spei$limits, breaks = y_info_spei$breaks) + 
          scale_color_manual(values = custom_colors, breaks = chrono_type_levels_ordered, name = "Wood Variable") +
          scale_shape_manual(values = custom_shapes, breaks = chrono_type_levels_ordered, name = "Wood Variable") +
          labs(
            title = paste("Correlation: Wood Variables and", current_spei_scale_label),
            subtitle = "SPEI calculated for each month (prev. year Oct-Dec, curr. year Jan-Dec).", # Adjusted subtitle
            x = "Month",
            y = "Correlation Coefficient (r)"
          ) +
          theme_bw() +
          theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=7),
            axis.text.y = element_text(size=8),
            axis.title = element_text(size=10),
            plot.title = element_text(size=12, hjust = 0.5, face="bold"),
            plot.subtitle = element_text(size=9, hjust = 0.5),
            legend.position = "bottom",
            legend.title = element_text(size=9, face="bold"),
            legend.text = element_text(size=8),
            strip.text.x = element_text(size = 9, face="bold"),
            strip.text.y = element_text(size = 9, face="bold"),
            panel.spacing = unit(0.8, "lines")
          ) +
          guides(color = guide_legend(title="Wood Variable"), shape = guide_legend(title="Wood Variable"))
        
        pdf_filename_spei <- paste0("Fig_Correlations_", current_spei_scale_label, "_By_Age.pdf")
        output_pdf_path_spei <- file.path(output_dir_spei_corr, pdf_filename_spei)
        
        ggsave(
          filename = output_pdf_path_spei,
          plot = spei_correlation_plot,
          width = 11,
          height = 9,
          units = "in",
          device = "pdf"
        )
        cat("    Plot for", current_spei_scale_label, "saved to:", output_pdf_path_spei, "\n")
      } # Fin del bucle for SPEI scales
    } # Fin del if (length(plot_spei_data_list) > 0)
  } # Fin del if (length(all_spei_data_by_scale) == 0)
} # Fin del if (length(spei_files_full_path) == 0)

cat("\n--- Proceso completado ---\n")
cat("\n--- Resumen de todas las advertencias de la sesión: ---\n")
if(!is.null(warnings()) && length(warnings()) > 0) print(warnings()) else cat("No hubo advertencias en esta sesión.\n")
