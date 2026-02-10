# ==============================================================================
# 01_climograma.R
# Autor: Marcos Marín-Martín
# Fecha: 2026-02-09
# Descripción: Generador de climograma comparando temperatura media mensual y
#              precipitación total.
# ==============================================================================

# --- 0. Cargar paquetes necesarios ---
# Ensure tidyverse and lubridate are installed
# install.packages(c("tidyverse", "lubridate"))
library(tidyverse)
library(lubridate)


# --- PARÁMETROS A MODIFICAR ---
# Rango temporal
YEAR_START <- 1985
YEAR_END   <- 2015

# Rutas de entrada
# Rutas de entrada
INPUT_TEMP_FILE   <- 'PLACEHOLDER/path/to/temp_data.txt'
INPUT_PRECIP_FILE <- 'PLACEHOLDER/path/to/precip_data.txt'

# Directorio de salida
OUTPUT_DIR <- 'PLACEHOLDER/path/to/output_dir'

# Configuración estética
COLOR_PRECIP <- "skyblue3"    # Azul para precip
COLOR_TEMP   <- "indianred2"  # Rojo para temp
# ------------------------------

# --- Define fixed font size in points and calculate geom_text size ---
fixed_font_size_pt <- 9
# ggplot's geom_text size parameter is roughly points / ggplot2::.pt
geom_text_size_value <- fixed_font_size_pt / ggplot2::.pt

# Function to process daily data to monthly
# (Function remains the same as previous version)
    group_by(month) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
  return(filtered_data)
}

process_climate_data <- function(path_data, data_type) {
  data <- read_tsv(path_data, col_types = cols()) %>%
    rename_with(~ "date", 1) %>%
    mutate(date = ymd(date)) %>%
    filter(!is.na(date)) %>%
    mutate(
      year = year(date),
      month = month(date)
    )
  value_col_name <- colnames(data)[2]
  if (is.na(value_col_name)) stop("Could not identify value column (expected second column).")
  filtered_data <- data %>%
    filter(year >= YEAR_START & year <= YEAR_END) %>%
    group_by(year, month) %>%
    summarise(
      value = if(data_type == "temp") mean(!!sym(value_col_name), na.rm = TRUE)
      else sum(!!sym(value_col_name), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(month) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
  return(filtered_data)
}

# Process temperature and precipitation data
# (Remains the same)
monthly_temp <- process_climate_data(INPUT_TEMP_FILE, "temp") %>%
  rename(temp = mean_value) %>% mutate(temp = round(temp, 1))
monthly_precip <- process_climate_data(INPUT_PRECIP_FILE, "precip") %>%
  rename(precip = mean_value) %>% mutate(precip = round(precip, 1))

# Combine data
# (Remains the same)
climograph_data <- left_join(monthly_precip, monthly_temp, by = "month") %>%
  mutate(month_name = factor(month.abb[month], levels = month.abb)) %>%
  tidyr::complete(month = 1:12, fill = list(precip = NA, temp = NA)) %>%
  mutate(month_name = factor(month.abb[month], levels = month.abb)) %>%
  filter(!is.na(month_name))

# --- Create climograph ---
# (ggplot object creation remains the same)
max_p <- max(climograph_data$precip, na.rm = TRUE)
y_limit_precip <- max_p * 1.15

climograph <- ggplot(climograph_data) +
  geom_col(aes(x = month_name, y = precip), fill = COLOR_PRECIP, alpha = 0.7, width = 0.7) +
  geom_line(aes(x = month_name, y = temp * 2, group = 1), color = COLOR_TEMP, linewidth = 1.2, lineend = "round") +
  geom_point(aes(x = month_name, y = temp * 2), color = COLOR_TEMP, size = 3) +
  scale_y_continuous(
    name = "Precipitation (mm)",
    sec.axis = sec_axis(~ . / 2, name = "Temperature (°C)"),
    limits = c(0, y_limit_precip),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Month", title = NULL, subtitle = NULL) + # Ensure titles are NULL
  theme_minimal() +
  theme(
    axis.title = element_text(size = fixed_font_size_pt),
    axis.text = element_text(size = fixed_font_size_pt),
    legend.text = element_text(size = fixed_font_size_pt),
    legend.title = element_text(size = fixed_font_size_pt),
    strip.text = element_text(size = fixed_font_size_pt),
    plot.tag = element_text(size = fixed_font_size_pt),
    plot.title = element_blank(), # Explicitly blank
    plot.subtitle = element_blank(), # Explicitly blank
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = y_limit_precip * 0.96, ymax = y_limit_precip, fill = COLOR_PRECIP, alpha = 0.7) +
  annotate("text", x = 1.7, y = y_limit_precip * 0.98, label = "Precipitation", hjust = 0, color = "black", size = geom_text_size_value) +
  annotate("segment", x = 0.5, xend = 1.5, y = y_limit_precip * 0.92, yend = y_limit_precip * 0.92, color = COLOR_TEMP, linewidth = 1.2) +
  annotate("text", x = 1.7, y = y_limit_precip * 0.92, label = "Temperature", hjust = 0, color = "black", size = geom_text_size_value)

# --- Save graph ---
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Define output file paths
pdf_output_path <- file.path(OUTPUT_DIR, "climograph_halfA4.pdf")
png_output_path <- file.path(OUTPUT_DIR, "climograph_halfA4.png")

# --- MODIFIED: Save PDF with half A4 width ---
half_a4_width_in <- 8.27 / 2
pdf_height_in <- 3.5 # Adjust this height if needed for better proportions

ggsave(filename = pdf_output_path,
       plot = climograph,
       device = "pdf",
       width = half_a4_width_in, # Set specific width
       height = pdf_height_in,     # Set adjusted height
       units = "in",
       useDingbats = FALSE)

# --- OPTIONAL: Adjust PNG dimensions proportionally ---
png_height_in <- pdf_height_in # Keep the same aspect ratio, adjust width to match
ggsave(filename = png_output_path,
       plot = climograph,
       device = "png",
       width = half_a4_width_in, # Match PDF width
       height = png_height_in,   # Match PDF height
       units = "in",
       dpi = 300)

# Show success message
message("Climograph successfully generated and saved to:\n", OUTPUT_DIR,
        "\nPDF dimensions (inches): ", round(half_a4_width_in, 2), " x ", pdf_height_in)