# ==============================================================================
# 04_heatmaps_correlaciones_por_edades_juntos_comprimido.R
# Autor: Marcos Marín-Martín
# Fecha: 2026-02-09
# Descripción: Genera un panel comparativo A4 transpuesto y comprimido heatmaps.
# ==============================================================================

# --- 0. Carga de Librerías ---
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)

# --- PARÁMETROS A MODIFICAR ---
OUTPUT_DIR          <- 'PLACEHOLDER/path/to/input_output_dir'
CORRELATION_FILE    <- "correlation_results.txt"

SIGNIFICANCE_LEVEL  <- 0.01
# ------------------------------




# --- 1. Definición de Rutas y Parámetros Globales ---
correlation_file_path <- file.path(OUTPUT_DIR, CORRELATION_FILE)
if (!file.exists(correlation_file_path)) stop(paste("El archivo de datos no se encontró:", correlation_file_path))

nonsignif_color <- "grey95"
color_neg_max <- "#083061"; color_zero <- "#FFFFFF"
positive_colors <- c("#ffffb2", "#ffeda0", "#FDE28A", "#FDE289", "#FDD774", "#FCC864", "#FBBA54", "#FBAE4A", "#FA9842", "#FA7836", "#F9582D", "#F64326", "#E31A1D", "#CF0E22", "#A50425", "#800226")
period_levels_monthly <- c(paste0("p", month.abb[10:12]), month.abb)
period_levels_extended <- c(period_levels_monthly, "DJF", "MAM", "JJA", "SON", "Annual")
chrono_type_levels_ordered <- c("TRW", "EW", "LW", "DBI", "EWBI", "LWBI")
climate_vars_classic <- c("tmax", "tmin", "precip")
spei_vars_ordered <- paste0("SPEI_", c(1, 3, 6, 12, 18, 24), "m")
all_climate_vars_ordered <- c(climate_vars_classic, spei_vars_ordered)
climate_var_titles <- c("Tmax", "Tmin", "Precip.", "SPEI-1", "SPEI-3", "SPEI-6", "SPEI-12", "SPEI-18", "SPEI-24")
age_group_titles <- c("Young", "Mid-Age", "Old")

x_labels_monthly_compact <- c("O*", "N*", "D*", "J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
x_labels_extended_compact <- c(x_labels_monthly_compact, "DJF", "MAM", "JJA", "SON", "Ann.")
names(x_labels_monthly_compact) <- period_levels_monthly
names(x_labels_extended_compact) <- period_levels_extended

# --- 2. Carga y Preparación de los Datos (sin cambios) ---
cat("--- Cargando y preparando los datos de correlación ---\n")
all_correlations_df <- read.table(correlation_file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
heatmap_data_prepared <- all_correlations_df %>%
  filter(period %in% period_levels_extended, climate_variable %in% all_climate_vars_ordered) %>%
  mutate(
    correlation_for_plot = ifelse(p_value < SIGNIFICANCE_LEVEL, correlation, NA_real_),
    age_group = factor(age_group, levels = c("yng", "mid", "old"), labels = age_group_titles),
    chrono_type = factor(chrono_type, levels = rev(chrono_type_levels_ordered)),
    period = factor(period, levels = period_levels_extended),
    climate_variable = factor(climate_variable, levels = all_climate_vars_ordered)
  ) %>%
  filter(!is.na(age_group) & !is.na(chrono_type) & !is.na(climate_variable) & !is.na(period))
cat("Datos preparados.\n")

# --- 3. CÁLCULO DE LA ESCALA DE COLOR GLOBAL (sin cambios) ---
cat("--- Calculando la escala de color global ---\n")
max_abs_corr <- max(abs(heatmap_data_prepared$correlation), na.rm = TRUE)
color_limit <- ceiling(max_abs_corr * 10) / 10
final_colors <- c(color_neg_max, color_zero, color_zero, positive_colors)
value_anchors_corr_scale <- c(-color_limit, -0.1, 0.1, seq(from = 0.1, to = color_limit, length.out = length(positive_colors)))
final_values_rescaled <- scales::rescale(value_anchors_corr_scale, from = c(-color_limit, color_limit))
all_possible_breaks <- seq(-1, 1, by = 0.2)
legend_breaks <- round(all_possible_breaks[all_possible_breaks >= -color_limit & all_possible_breaks <= color_limit], 1)

# --- 4. GENERACIÓN DE TODOS LOS 27 PLOTS EN MEMORIA ---
cat("--- Generando los 27 heatmaps individuales en memoria ---\n")
all_plots <- list()
font_base_size <- 7
common_theme <- theme_minimal(base_size = font_base_size) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank(),
        plot.margin = unit(c(1, 2, 1, 2), "pt"))

for (current_clim_var in all_climate_vars_ordered) {
  all_plots[[current_clim_var]] <- list()
  is_extended <- current_clim_var %in% climate_vars_classic
  for (current_age_group in levels(heatmap_data_prepared$age_group)) {
    data_for_plot <- heatmap_data_prepared %>% 
      filter(age_group == current_age_group, climate_variable == current_clim_var)
    if (!is_extended) {
      data_for_plot <- data_for_plot %>% filter(period %in% period_levels_monthly)
    }
    p <- ggplot(data_for_plot, aes(x = period, y = chrono_type, fill = correlation_for_plot)) +
      geom_tile(color = "grey80", linewidth = 0.15) +
      scale_fill_gradientn(colors = final_colors, values = final_values_rescaled, limit = c(-color_limit, color_limit), name = "r", na.value = nonsignif_color, breaks = legend_breaks) +
      coord_fixed(ratio = 1) + labs(x = NULL, y = NULL) + common_theme
    if (is_extended) {
      p <- p + scale_x_discrete(drop = FALSE, labels = x_labels_extended_compact) +
        geom_vline(xintercept = length(period_levels_monthly) + 0.5, color = "gray50", linetype = "dashed", linewidth=0.2)
    } else {
      p <- p + scale_x_discrete(drop = FALSE, labels = x_labels_monthly_compact)
    }
    all_plots[[current_clim_var]][[current_age_group]] <- p
  }
}
cat("Todos los heatmaps generados.\n")

# --- 5. ENSAMBLAJE DEL PANEL TRANSPUESTO (3x9) ---
cat("--- Ensamblando el panel comparativo transpuesto ---\n")
plots_for_grid <- list()
for (clim_var in all_climate_vars_ordered) {
  for (age_group in levels(heatmap_data_prepared$age_group)) {
    plots_for_grid[[length(plots_for_grid) + 1]] <- all_plots[[clim_var]][[age_group]]
  }
}
for (i in seq_along(plots_for_grid)) {
  row_idx <- floor((i - 1) / 3) + 1; col_idx <- ((i - 1) %% 3) + 1
  if (row_idx < 9) plots_for_grid[[i]] <- plots_for_grid[[i]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  if (col_idx > 1) plots_for_grid[[i]] <- plots_for_grid[[i]] + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  if (row_idx == 1) {
    plots_for_grid[[i]] <- plots_for_grid[[i]] + labs(title = age_group_titles[col_idx]) +
      theme(plot.title = element_text(size = font_base_size + 2, face = "bold", hjust = 0.5))
  }
  if (col_idx == 1) {
    plots_for_grid[[i]] <- plots_for_grid[[i]] + labs(y = climate_var_titles[row_idx]) +
      theme(axis.title.y = element_text(angle=0, vjust=0.5, face="bold", size=font_base_size+1))
  }
}

final_plot <- wrap_plots(plots_for_grid, ncol = 3, guides = 'collect') & 
  theme(
    legend.position = 'bottom',
    # *** CAMBIO CLAVE: Ancho de la leyenda reducido para que quepa ***
    legend.key.width = unit(2.8, "cm"), 
    legend.key.height = unit(0.4, "cm")
  )

# --- 6. GUARDAR EL PDF FINAL ---
output_filename <- file.path(OUTPUT_DIR, "Heatmaps_Correlations_AllAges_A4_Transposed_Final_Polished.pdf")
cat("--- Guardando el panel comparativo final en:", output_filename, "---\n")
ggsave(
  filename = output_filename,
  plot = final_plot,
  width = 210, 
  height = 297,
  units = "mm",
  device = "pdf",
  limitsize = FALSE
)
cat("\n--- Proceso completado ---\n")