# ==============================================================================
# 09_plot_sss_trwi_depth.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Standard Chronology Quality Visualisation (SSS + RWI + Depth).
#   Generates a publication-quality 3-panel plot to assess the reliability of a 
#   standardised tree-ring chronology.
#
#   Panels:
#   1. Subsample Signal Strength (SSS): Monitors the decline in signal reliability 
#      back in time as sample size decreases. Adds a vertical line where SSS drops 
#      below a critical threshold (e.g., 0.85).
#   2. Tree-Ring Width Index (RWI): The standardised chronology itself (Residual or Standard).
#   3. Sample Depth: Number of tree series contributing to the chronology over time.
#
#   Inputs:
#   - RWL file (for SSS calculation).
#   - Chronology TXT file (for RWI and Depth).
#
#   Outputs:
#   - PDF Plot.
# ==============================================================================

library(dplR)


# --- PARÁMETROS A MODIFICAR ---
# Estética
COL_SSS         <- "orangered4"  # Color para SSS
COL_TRW         <- "steelblue1"  # Color para TRW Index
COL_SAMP_DEPTH  <- "orangered4"  # Color para Sample Depth
LWD_SSS         <- 1             # Ancho de línea para SSS
LWD_TRW         <- 2             # Ancho de línea para TRW Index
LWD_SAMP_DEPTH  <- 1             # Ancho de línea para Sample Depth

# Salida
OUTPUT_FILE     <- "plot_sss_trw_startyear_vline_fix.pdf"
OUTPUT_PATH     <- "PLACEHOLDER/path/to/output_dir/"

# Dimensiones y Configuración
PLOT_WIDTH      <- 11         # Ancho del gráfico
PLOT_HEIGHT     <- 6          # MANTENEMOS LA ALTURA REDUCIDA
FONT_SIZE       <- 10         # Tamaño de la fuente
SSS_THRESHOLD   <- 0.85       # Umbral SSS para la LÍNEA VERTICAL
START_PLOT_YEAR <- 1443       # Año de inicio para mostrar en el gráfico

# Línea vertical SSS
VLINE_COL       <- "grey20"    # Color de la línea vertical SSS
VLINE_LTY       <- 2           # Tipo de línea vertical SSS (discontinua)
VLINE_LWD       <- 1.5         # Ancho de la línea vertical SSS

# Rutas de entrada
INPUT_RWL_FILE  <- 'PLACEHOLDER/path/to/input_data.rwl'
INPUT_TXT_FILE  <- 'PLACEHOLDER/path/to/chronology.txt'
# ------------------------------


# ----------------- CARGA DE DATOS ----------------- #

# Cargar datos de la cronología .rwl
data_rwl <- read.rwl(INPUT_RWL_FILE)

# Calcular SSS
sss_results <- sss(data_rwl)
# --- CORREGIDO: Usar names() para extraer años del vector SSS ---
years_sss <- as.numeric(names(sss_results)) 

# --- MODIFICADO: Encontrar el año umbral pero NO filtrar por él ---
# Identificar años donde SSS >= umbral (solo para encontrar el primero)
# Añadido !is.na() para seguridad
years_above_threshold <- years_sss[!is.na(sss_results) & sss_results >= SSS_THRESHOLD]
if (length(years_above_threshold) > 0) {
  first_year_threshold <- min(years_above_threshold)
  print(paste("El umbral SSS de", SSS_THRESHOLD, "se alcanza por primera vez en:", first_year_threshold))
} else {
  first_year_threshold <- NA # O manejar de otra forma si nunca se alcanza
  # La advertencia ya se produjo según tu output, indicando que esto es correcto
  warning(paste("El umbral SSS de", SSS_THRESHOLD, "nunca se alcanza o no hay datos SSS válidos."))
}

# Filtrar datos SSS para el gráfico (desde START_PLOT_YEAR)
sss_plot_years <- years_sss[years_sss >= START_PLOT_YEAR]
# --- CORREGIDO: Subsetting correcto para vector con nombres ---
sss_plot_values <- sss_results[as.character(sss_plot_years)]


# Cargar datos del archivo .txt
data_txt <- read.table(INPUT_TXT_FILE, header = TRUE)

# --- MODIFICADO: Filtrar datos TXT solo por START_PLOT_YEAR ---
filtered_data_txt <- data_txt[data_txt$year >= START_PLOT_YEAR, ]
if (nrow(filtered_data_txt) == 0) {
  stop(paste("No hay datos en el archivo .txt a partir del año", START_PLOT_YEAR))
}
years_plot <- filtered_data_txt$year
samp_depth_plot <- filtered_data_txt$samp.depth
trw_index_plot <- filtered_data_txt$res

# --- MODIFICADO: Rango X común para todos los paneles ---
x_range_plot <- c(START_PLOT_YEAR, max(years_plot))

# Encontrar el primer año redondo (50 o 00) DENTRO DEL RANGO DEL GRÁFICO
first_round_year_plot <- ceiling(START_PLOT_YEAR / 50) * 50

# ----------------- GENERACIÓN DEL GRÁFICO ----------------- #
pdf(file.path(OUTPUT_PATH, OUTPUT_FILE), width = PLOT_WIDTH, height = PLOT_HEIGHT)
layout(matrix(c(1, 2, 3), nrow = 3), heights = c(3, 4, 3))

# Configuración común: sin bordes automáticos y con márgenes exteriores
par(xaxs = "i", yaxs = "i", bty = "n", cex = FONT_SIZE / 10, oma = c(1, 1, 1, 1))

# --- Panel SSS (sin bordes) ---
par(mar = c(0, 4, 0.5, 4))
# Asegurar que ylim tiene sentido incluso si sss_plot_values está vacío o todo NA
min_y_sss <- min(c(sss_plot_values, SSS_THRESHOLD), na.rm = TRUE)
max_y_sss <- 1 # SSS no puede ser > 1
if (!is.finite(min_y_sss)) min_y_sss <- SSS_THRESHOLD - 0.05 # Valor por defecto si no hay datos
plot(sss_plot_years, sss_plot_values, type = "l", col = COL_SSS, lwd = LWD_SSS,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n",
     xlim = x_range_plot,
     ylim = c(min_y_sss, max_y_sss)) # Rango Y ajustado
axis(4, at = pretty(c(min_y_sss, max_y_sss)), lwd = 0, lwd.ticks = 1, cex.axis = FONT_SIZE / 10) # Eje automático
mtext("SSS", side = 4, line = 2, cex = FONT_SIZE / 10)
# --- Añadir línea vertical ---
if (!is.na(first_year_threshold) && first_year_threshold >= START_PLOT_YEAR) {
  abline(v = first_year_threshold, col = VLINE_COL, lty = VLINE_LTY, lwd = VLINE_LWD)
}

# --- Panel TRW (sin bordes) ---
par(mar = c(0, 4, 0, 4))
plot(years_plot, trw_index_plot, type = "l", col = COL_TRW, lwd = LWD_TRW,
     xaxt = "n", xlab = "", ylab = "", bty = "n",
     xlim = x_range_plot)
axis(2, lwd = 0, lwd.ticks = 1, cex.axis = FONT_SIZE / 10)
mtext("TRW Index", side = 2, line = 2, cex = FONT_SIZE / 10)
abline(h = mean(trw_index_plot, na.rm=TRUE), col = "black", lty = 2, lwd = 1)
# --- Añadir línea vertical ---
if (!is.na(first_year_threshold) && first_year_threshold >= START_PLOT_YEAR) {
  abline(v = first_year_threshold, col = VLINE_COL, lty = VLINE_LTY, lwd = VLINE_LWD)
}

# --- Panel Sample Depth (sin bordes) ---
par(mar = c(4, 4, 0, 4))
plot(years_plot, samp_depth_plot, type = "l", col = COL_SAMP_DEPTH, lwd = LWD_SAMP_DEPTH,
     xlab = "Years", yaxt = "n", ylab = "", bty = "n",
     xlim = x_range_plot,
     ylim = c(0, max(samp_depth_plot, na.rm=TRUE)))
axis(4, at = pretty(c(0, max(samp_depth_plot, na.rm=TRUE))), lwd = 0, lwd.ticks = 1, cex.axis = FONT_SIZE / 10)
axis(1, at = seq(first_round_year_plot, max(years_plot), by = 50), lwd = 0, lwd.ticks = 1, cex.axis = FONT_SIZE / 10)
mtext("Sample Depth", side = 4, line = 2, cex = FONT_SIZE / 10)
# --- Añadir línea vertical ---
if (!is.na(first_year_threshold) && first_year_threshold >= START_PLOT_YEAR) {
  abline(v = first_year_threshold, col = VLINE_COL, lty = VLINE_LTY, lwd = VLINE_LWD)
}

# Añadir el rectángulo negro alrededor de los tres paneles
par(fig = c(0, 1, 0, 1), new = TRUE)
plot(0, 0, type = "n", xlab = "", ylab = "", axes = FALSE)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], border = "black", lwd = 1)

dev.off()

print(paste("Gráfico corregido guardado en:", file.path(OUTPUT_PATH, OUTPUT_FILE)))
