# ==============================================================================
# 05_plot_reconst_txt.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Publication-Ready Reconstruction Plotter.
#   A highly polished plotting utility designed to generate final figures for 
#   publication from a simple text file input.
#
#   Aesthetics & Features:
#   - Clean, minimal design with 'Rubik' font support (if available) or defaults.
#   - Visualizes:
#     - Reconstructed series (Cyan).
#     - Observed instrumental series (Red Dashed).
#     - Low-frequency signal (11-year Rolling Mean, Grey).
#     - Uncertainty Ribbon based on RMSE (Purple).
#   - Extreme Events: Automatically identifies and labels the wettest and driest 
#     years per century or for the whole period.
#
#   Inputs:
#   - Reconstruction Result Text File.
#
#   Outputs:
#   - High-Quality PNG Figure (600 dpi).
# ==============================================================================

# Cargar las librerías necesarias
library(ggplot2)
library(zoo)
library(extrafont)
library(dplyr)


# --- PARÁMETROS A MODIFICAR ---
INPUT_FILE           <- "PLACEHOLDER/path/to/reconstruction.txt"
OUTPUT_FILE_PNG      <- "PLACEHOLDER/path/to/precipitacion_plot_v2.png"
ROLLING_WINDOW_YEARS <- 11
# ------------------------------

# Cargar la fuente Rubik (debe estar instalada previamente en el sistema)
loadfonts(device = "win", quiet = TRUE)

# Cargar los datos desde el archivo
data <- read.table(INPUT_FILE, header = TRUE)

# Renombrar las columnas para mayor claridad
names(data) <- c("Year", "std", "TRW", "samp.depth", "Precipitación_real", "Precipitación_reconstruida", "media_movil_11")

# Calcular la media móvil centrada de 11 años
data$media_movil_11 <- zoo::rollmean(data$Precipitación_reconstruida, ROLLING_WINDOW_YEARS, fill = NA, align = "center")

# Calcular métricas de error
calibration_mean <- mean(data$Precipitación_reconstruida, na.rm = TRUE)
rmse <- sqrt(mean((data$Precipitación_reconstruida - data$Precipitación_real)^2, na.rm = TRUE))

# Identificar extremos por siglo
data <- data %>%
  mutate(Century = floor(Year/100)*100) %>%
  group_by(Century) %>%
  mutate(
    min_century = ifelse(Precipitación_reconstruida == min(Precipitación_reconstruida, na.rm = TRUE), Year, NA),
    max_century = ifelse(Precipitación_reconstruida == max(Precipitación_reconstruida, na.rm = TRUE), Year, NA)
  ) %>%
  ungroup()

# Identificar los x años más secos y los x años más húmedos de todo el período
driest_years <- data %>% arrange(Precipitación_reconstruida) %>% head(8)
wettest_years <- data %>% arrange(desc(Precipitación_reconstruida)) %>% head(7)

# Crear etiquetas para la leyenda
leyenda_rmse <- paste0("RMSE = ", round(rmse, 1), " mm")

# Crear el gráfico con todas las modificaciones
ggplot(data, aes(x = Year)) +
  geom_line(aes(y = Precipitación_reconstruida, color = "Reconstructed"), size = 0.3) +
  geom_line(aes(y = Precipitación_real, color = "Observed"), linetype = "dashed", size = 0.5) +
  geom_line(aes(y = media_movil_11, color = "11 years rolling mean"), size = 0.8) +
  geom_ribbon(aes(ymin = Precipitación_reconstruida - rmse, 
                  ymax = Precipitación_reconstruida + rmse, 
                  fill = leyenda_rmse), alpha = 0.2) +
  geom_point(data = filter(data, !is.na(min_century)), 
             aes(y = Precipitación_reconstruida), color = "red", size = 1.5) +
  geom_point(data = filter(data, !is.na(max_century)), 
             aes(y = Precipitación_reconstruida), color = "blue", size = 1.5) +
  geom_point(data = driest_years, aes(y = Precipitación_reconstruida), color = "red", size = 1.5) +
  geom_point(data = wettest_years, aes(y = Precipitación_reconstruida), color = "blue", size = 1.5) +
  geom_text(data = filter(data, !is.na(min_century)), 
            aes(y = Precipitación_reconstruida, label = min_century),
            family = "Rubik", angle = 45, hjust = 1.1, vjust = 1.5, size = 2.5, color = "black") +
  geom_text(data = filter(data, !is.na(max_century)), 
            aes(y = Precipitación_reconstruida, label = max_century),
            family = "Rubik", angle = 45, hjust = -0.1, vjust = -0.5, size = 2.5, color = "black") +
  geom_text(data = driest_years, aes(y = Precipitación_reconstruida, label = Year),
            family = "Rubik", angle = 45, hjust = 1.1, vjust = 1.5, size = 2.5, color = "black") +
  geom_text(data = wettest_years, aes(y = Precipitación_reconstruida, label = Year),
            family = "Rubik", angle = 45, hjust = -0.1, vjust = -0.5, size = 2.5, color = "black") +
  scale_color_manual(values = c(
    "Observed" = "#E41A1C",
    "Reconstructed" = "cyan3",
    "11 years rolling mean" = "gray50")) +
  scale_fill_manual(values = c("#984EA3")) +
  labs(
    title = "Precipitations (13 August previous year - 1 July current year)",
    x = "Year",
    y = "Precipitation (mm)",
    color = "",
    fill = "Uncertainty") +
  theme_minimal(base_family = "Rubik") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank()) +
  guides(
    color = guide_legend(order = 1, nrow = 1),
    fill = guide_legend(order = 2)) +
  scale_x_continuous(breaks = seq(1500, 2020, by = 50))

# Guardar el gráfico con calidad alta
ggsave(OUTPUT_FILE_PNG, 
       width = 30, 
       height = 15, 
       units = "cm", 
       dpi = 600)