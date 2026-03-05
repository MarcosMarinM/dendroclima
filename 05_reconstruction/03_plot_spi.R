# ==============================================================================
# 04_plot_spi.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   SPI Reconstruction Visualizer.
#   Specialised plotting script for Standardised Precipitation Index (SPI) reconstructions.
#   Visualizes the temporal evolution of drought/wetness with specific SPI thresholds.
#
#   Visualisation Features:
#   - Reconstruction line with uncertainty ribbon (RMSE).
#   - Observed SPI comparison (dashed line).
#   - 11-Year Moving Average.
#   - Critical Threshold Lines:
#     - Horizontal lines for SPI > 1 (Wet) and SPI < -1 (Dry).
#     - Annotated labels for these thresholds.
#
#   Inputs:
#   - SPI Reconstruction Data File.
#
#   Outputs:
#   - High-resolution Plot (.png, 600dpi).
#   - Vector Graphics Plot (.pdf).
# ==============================================================================


library(ggplot2)
library(zoo)
library(dplyr)

# --- PARÁMETROS A MODIFICAR ---
INPUT_FILE           <- "PLACEHOLDER/path/to/spi_reconstruction.txt"
OUTPUT_FILE_PNG      <- "PLACEHOLDER/path/to/spi_plot.png"
OUTPUT_FILE_PDF      <- "PLACEHOLDER/path/to/spi_plot.pdf"
ROLLING_WINDOW_YEARS <- 11  # Media móvil (en años)
SPI_THRESHOLDS       <- c(-2, -1.5, -1, 1, 1.5, 2)  # Umbrales para clasificar SPI
# ------------------------------

# Cargar los datos
data <- read.table(INPUT_FILE, header = TRUE)

# Filtrar datos a partir del primer año con valores en 'reconstruccion' o 'spi_320d_0701'
data <- data %>% filter(!is.na(reconstruccion) | !is.na(spi_320d_0701))

# Calcular la media móvil centrada
data$media_movil_11 <- zoo::rollmean(data$reconstruccion, ROLLING_WINDOW_YEARS, fill = NA, align = "center")

# Calcular RMSE entre reconstrucción y observación
rmse <- sqrt(mean((data$reconstruccion - data$spi_320d_0701)^2, na.rm = TRUE))

# Crear el gráfico
p <- ggplot(data, aes(x = year)) +
  
  # Margen de error con RMSE
  geom_ribbon(aes(ymin = reconstruccion - rmse, ymax = reconstruccion + rmse), fill = "lightskyblue3", alpha = 0.2) +
  
  # Líneas de SPI reconstruido, observado y media móvil
  geom_line(aes(y = reconstruccion, color = "Reconstrucción"), size = 0.8) +
  geom_line(aes(y = spi_320d_0701, color = "Observación"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = media_movil_11, color = "Media móvil 11 años"), size = 0.8) +
  
  # Líneas horizontales en 1 y -1 con etiquetas
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue", size = 0.3) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "red", size = 0.3) +
  annotate("text", x = max(data$year, na.rm = TRUE) + 10, y = 1.2, label = "Años húmedos (SPI > 1)", color = "blue", hjust = 1, size = 4) +
  annotate("text", x = max(data$year, na.rm = TRUE) + 10, y = -1.2, label = "Años secos (SPI < -1)", color = "red", hjust = 1, size = 4) +
  
  # Configuración de colores
  scale_color_manual(values = c("Observación" = "indianred2", "Reconstrucción" = "steelblue2", "Media móvil 11 años" = "gray50")) +
  
  # Títulos y etiquetas
  labs(
    title = "Reconstrucción del SPI (16 agosto año previo - 30 junio año actual)",
    x = "Año",
    y = "SPI",
    color = "") +
  
  # Estilos de la gráfica
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  
  # Eje X más ajustado
  scale_x_continuous(limits = c(min(data$year, na.rm = TRUE), max(data$year, na.rm = TRUE)), 
                     breaks = seq(min(data$year, na.rm = TRUE), max(data$year, na.rm = TRUE), by = 50))

# Guardar el gráfico en PNG
ggsave(OUTPUT_FILE_PNG, plot = p, width = 30, height = 15, units = "cm", dpi = 600)

# Guardar el gráfico en PDF
ggsave(OUTPUT_FILE_PDF, plot = p, width = 30, height = 15, units = "cm")

# Mostrar el gráfico en R
print(p)
