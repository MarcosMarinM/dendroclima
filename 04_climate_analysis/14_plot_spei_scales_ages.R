# ==============================================================================
# 04_grafico_escalas_spei_edades.R
# Autor: Marcos Marín-Martín
# Fecha: 2026-02-09
# Descripción: Calcula y visualiza correlaciones entre cronologías y SPEI.
#              Análisis por grupos de edad y escalas temporales.
# ==============================================================================

# --- 0. Configuración y Carga de Librerías ---
cat("--- 0. Cargando librerías necesarias ---\n")
# --- 0. Configuración y Carga de Librerías ---
cat("--- 0. Cargando librerías necesarias ---\n")
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(here)


# --- PARÁMETROS A MODIFICAR ---
# Rutas
# --- PARÁMETROS A MODIFICAR ---
# Rutas
CHRONO_DIR <- "PLACEHOLDER/path/to/chrono_dir"
SPEI_DIR   <- "PLACEHOLDER/path/to/spei_dir"
OUTPUT_DIR <- "PLACEHOLDER/path/to/output_dir"

# Parámetros de Análisis
MIN_OVERLAP_YEARS    <- 30
SIGNIFICANCE_LEVEL   <- 0.05
MONTH_TO_PLOT        <- "Jun"
CHRONO_TYPES_TO_PLOT <- c("TRW", "DBI")

# Estética
AGE_GROUP_COLORS <- c("Young" = "forestgreen", "Mid-Age" = "dodgerblue", "Old" = "firebrick")
AGE_GROUP_SHAPES <- c("Young" = 16, "Mid-Age" = 17, "Old" = 15)
AGE_GROUP_LEVELS <- c("Young", "Mid-Age", "Old")
# ------------------------------

# --- 2. Carga y Procesamiento de Cronologías ---
chrono_files <- list.files(CHRONO_DIR, pattern = "\\.txt$", full.names = TRUE)

chrono_list <- lapply(chrono_files, function(file_path) {
  # Extraer metadatos del nombre del archivo
  file_name <- basename(file_path)
  # Usamos regex para capturar tipo y edad de forma robusta
  meta <- stringr::str_match(file_name, "IBE_(\\w+)_(\\w+)\\.txt")
  
  if (is.na(meta[1,1])) return(NULL) # Si el nombre no coincide, lo saltamos
  
  chrono_type <- meta[1, 2]
  age_group_short <- meta[1, 3]
  
  # Estandarizar nombres de grupos de edad
  age_group <- case_when(
    age_group_short == "yng" ~ "Young",
    age_group_short == "mid" ~ "Mid-Age",
    age_group_short == "old" ~ "Old",
    TRUE ~ age_group_short
  )
  
  # Leer los datos, seleccionando solo las columnas necesarias
  df <- readr::read_table(file_path, col_types = readr::cols()) %>%
    select(year, res)
  
  # Devolver una lista con los datos y metadatos
  list(
    type = chrono_type,
    age = age_group,
    data = df
  )
})
# Eliminar cualquier archivo que no coincidiera con el patrón
chrono_list <- purrr::compact(chrono_list)
cat("Cargadas", length(chrono_list), "cronologías.\n")


# --- 3. Carga y Transformación de Datos SPEI ---
cat("--- 3. Cargando y transformando datos SPEI a formato 'tidy' ---\n")
spei_files <- list.files(SPEI_DIR, pattern = "\\.txt$", full.names = TRUE)
month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

spei_tidy_df <- purrr::map_dfr(spei_files, function(file_path) {
  # Extraer la escala del nombre del archivo
  scale <- as.numeric(stringr::str_extract(basename(file_path), "\\d+"))
  
  # Leer el archivo, nombrando las columnas de mes V1, V2, etc.
  df_wide <- readr::read_table(file_path, col_names = FALSE, col_types = readr::cols())
  names(df_wide) <- c("year", month_names)
  
  # Transformar de formato ancho a largo (tidy)
  df_long <- df_wide %>%
    tidyr::pivot_longer(
      cols = all_of(month_names),
      names_to = "period",
      values_to = "spei_value"
    ) %>%
    mutate(spei_scale = scale)
  
  return(df_long)
})
cat("Datos SPEI transformados. Total de filas:", nrow(spei_tidy_df), "\n")


# --- 4. Cálculo de Correlaciones ---
cat("--- 4. Calculando correlaciones (esto puede tardar un poco)... ---\n")
results_list <- list()
total_calcs <- length(chrono_list) * length(unique(spei_tidy_df$spei_scale)) * 12
progress_count <- 0

# Iteramos sobre cada cronología cargada
for (chrono in chrono_list) {
  chrono_df <- chrono$data
  
  # Iteramos sobre cada combinación única de escala y mes del SPEI
  # Usar split() es más eficiente que un doble bucle for
  spei_groups <- split(spei_tidy_df, ~ spei_scale + period)
  
  for (spei_subset in spei_groups) {
    progress_count <- progress_count + 1
    cat(sprintf("  Calculando %d de %d...\r", progress_count, total_calcs))
    
    # Unir por año para encontrar el período de solapamiento
    merged_data <- dplyr::inner_join(chrono_df, spei_subset, by = "year")
    
    # Comprobación de seguridad: ¿hay suficientes datos?
    if (nrow(merged_data) < MIN_OVERLAP_YEARS) {
      next # Saltar al siguiente si no hay suficientes años
    }
    
    # Calcular correlación de Pearson
    cor_result <- cor.test(~ res + spei_value, data = merged_data, method = "pearson")
    
    # Almacenar el resultado
    results_list[[length(results_list) + 1]] <- tibble::tibble(
      chrono_type = chrono$type,
      age_group = chrono$age,
      spei_scale = unique(spei_subset$spei_scale),
      period = unique(spei_subset$period),
      correlation = cor_result$estimate,
      p_value = cor_result$p.value,
      n_overlap = nrow(merged_data)
    )
  }
}
cat("\nCálculo de correlaciones completado.\n")

# Combinar todos los resultados en un solo data.frame
correlations_df <- dplyr::bind_rows(results_list)


# --- 5. Guardado de Resultados ---
cat("--- 5. Guardando tabla de correlaciones ---\n")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
output_rds_path <- file.path(OUTPUT_DIR, "chronology_spei_correlations.rds")
saveRDS(correlations_df, file = output_rds_path)
cat("Resultados guardados en:", output_rds_path, "\n")


# --- 6. Visualización ---
cat("--- 6. Generando el gráfico de respuesta al SPEI ---\n")

# Preparar los datos para el plot (filtrar y añadir columnas de ayuda)
plot_data <- correlations_df %>%
  filter(
    period == MONTH_TO_PLOT,
    chrono_type %in% CHRONO_TYPES_TO_PLOT
  ) %>%
  mutate(
    is_significant = p_value < SIGNIFICANCE_LEVEL,
    age_group = factor(age_group, levels = AGE_GROUP_LEVELS),
    period = factor(period, levels = month_names) # Para ordenar meses si graficas varios
  ) %>%
  arrange(chrono_type, age_group, spei_scale)

# Preparar datos para los segmentos de línea significativos
line_segments_df <- plot_data %>%
  group_by(chrono_type, age_group) %>%
  mutate(
    x_end = lead(spei_scale),
    y_end = lead(correlation),
    next_is_significant = lead(is_significant)
  ) %>%
  ungroup() %>%
  filter(is_significant & next_is_significant)

# Crear el objeto ggplot
final_plot <- ggplot(plot_data, aes(x = spei_scale, y = correlation, color = age_group, group = age_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.5, alpha = 0.4) +
  geom_segment(data = line_segments_df, aes(xend = x_end, yend = y_end), linewidth = 1.2) +
  geom_point(aes(alpha = is_significant, size = is_significant, shape = age_group)) +
  facet_wrap(~ chrono_type, ncol = 1, scales = "free_y") +
  scale_y_continuous(name = "Correlation Coefficient (r)") +
  scale_x_continuous(name = "SPEI Timescale (months)", breaks = seq(0, 24, by = 3)) +
  scale_color_manual(name = "Age Group", values = AGE_GROUP_COLORS) +
  scale_shape_manual(name = "Age Group", values = AGE_GROUP_SHAPES) +
  scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.3), guide = "none") +
  scale_size_manual(values = c("TRUE" = 3.0, "FALSE" = 1.5), guide = "none") +
  labs(
    title = paste("Age-dependent Response to", MONTH_TO_PLOT, "SPEI Timescale"),
    subtitle = paste("Highlighted correlations are significant (p <", SIGNIFICANCE_LEVEL, ")"),
    color = "Age Group",
    shape = "Age Group"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))

# Guardar el gráfico
output_plot_path <- file.path(OUTPUT_DIR, paste0("SPEI_Response_Plot_", MONTH_TO_PLOT, ".pdf"))
ggsave(output_plot_path, plot = final_plot, width = 8, height = 7, device = "pdf")

cat("--- ¡Proceso completado! ---\n")
cat("Gráfico guardado en:", output_plot_path, "\n")