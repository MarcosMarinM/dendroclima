# ==============================================================================
# 09_plot_sss_trwi_rbar_depth.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Comprehensive Chronology Diagnostics (EPS/SSS + Rbar + RWI + Depth).
#   An advanced visualisation script creates a composite plot for detailed 
#   chronology assessment. It combines running signal statistics with the main 
#   chronology series.
#
#   Key Metrics Visualised:
#   - SSS (Subsample Signal Strength): Reliability over time.
#   - Rbar (Running Mean Correlation): Internal coherency of the series windowed over time.
#   - TRW Index: The standardised residual chronology.
#   - Sample Depth: Sample replication.
#
#   Inputs:
#   - RWI file.
#   - Chronology TXT file.
#   - Reconstruction/Analysis output (optional, for specific SSS data).
#
#   Outputs:
#   - "Final Figure" style PDF for publication or deep analysis.
# ==============================================================================

library(dplR)

# ----------------- PARÁMETROS ----------------- #

# --- PARÁMETROS A MODIFICAR ---
# Estética
COL_SSS         <- "coral2"
COL_RBAR        <- "gold3"
COL_TRW         <- "steelblue1"
COL_SAMP_DEPTH  <- "aquamarine4"

LWD_SSS         <- 1.5
LWD_RBAR        <- 1.5
LWD_TRW         <- 1.5
LWD_SAMP_DEPTH  <- 1.5

# Salida y Configuración
OUTPUT_FILE     <- "Fig2_FINAL_FINAL.pdf"
OUTPUT_PATH     <- "PLACEHOLDER/path/to/output_dir/"
PLOT_WIDTH      <- 11
PLOT_HEIGHT     <- 7
FONT_SIZE       <- 10
SSS_THRESHOLD   <- 0.85
START_PLOT_YEAR <- 1443
WINDOW_RBAR     <- 51

# Rutas de entrada
INPUT_RWI_FILE            <- 'PLACEHOLDER/path/to/input_data.rwi'
INPUT_TXT_FILE            <- 'PLACEHOLDER/path/to/chronology.txt'
INPUT_RECONSTRUCTION_FILE <- 'PLACEHOLDER/path/to/reconstruction_data.txt'
# ------------------------------


# ----------------- CARGA DE DATOS ----------------- #
data_rwi <- read.rwl(INPUT_RWI_FILE)
data_txt <- read.table(INPUT_TXT_FILE, header = TRUE)
data_reconstruction <- read.table(INPUT_RECONSTRUCTION_FILE, header = TRUE)

# --- CÁLCULO DE MÉTRICAS Y PREPARACIÓN --- #
running_stats <- rwi.stats.running(data_rwi, window.length = WINDOW_RBAR, window.overlap = WINDOW_RBAR - 1)

sss_above_threshold_df <- subset(data_reconstruction, sss >= SSS_THRESHOLD & !is.na(sss))
first_year_threshold <- if (nrow(sss_above_threshold_df) > 0) min(sss_above_threshold_df$year) else NA

x_range_plot <- c(START_PLOT_YEAR, max(data_txt$year))
sss_to_plot <- subset(data_reconstruction, year >= START_PLOT_YEAR) 
rbar_to_plot <- subset(running_stats, mid.year >= START_PLOT_YEAR)
txt_to_plot <- subset(data_txt, year >= START_PLOT_YEAR)

# ----------------- GENERACIÓN DEL GRÁFICO ----------------- #
pdf(file.path(OUTPUT_PATH, OUTPUT_FILE), width = PLOT_WIDTH, height = PLOT_HEIGHT)

### MODIFICADO: Layout de 2 paneles con proporción 60% / 40% ###
layout(matrix(1:2, nrow = 2), heights = c(3, 2))
par(xaxs = "i", yaxs = "i", bty = "n", cex = FONT_SIZE / 10, oma = c(1, 1, 1, 1))

# --- Panel 1: TRW Index (arriba, grande) ---
par(mar = c(0, 4.5, 1.5, 4.5))
plot(txt_to_plot$year, txt_to_plot$res, type = "l", col = COL_TRW, xaxt = "n", xlab = "", ylab = "", bty = "n", xlim = x_range_plot, lwd=LWD_TRW)
axis(2, lwd = 0, lwd.ticks = 1); mtext("TRW Index (res)", side = 2, line = 3)
abline(h = 0, col = "black", lty = 1, lwd = 0.5)
### MODIFICADO: Estilo de la línea vertical ###
if (!is.na(first_year_threshold)) { abline(v = first_year_threshold, col = COL_SSS, lty = "dotted", lwd = 1.5) }

# --- Panel 2: Métricas de Calidad (abajo, compacto) ---
par(mar = c(4, 4.5, 0.5, 4.5)) 

# Dibujar SSS y Rbar en el eje izquierdo
plot(sss_to_plot$year, sss_to_plot$sss, type = "l", col = COL_SSS, lwd = LWD_SSS, 
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", 
     xlim = x_range_plot, ylim = c(0, 1))
### MODIFICADO: Eje con saltos de 0.2 ###
axis(2, at = seq(0, 1, by = 0.2), lwd = 0, lwd.ticks = 1); mtext("SSS / Running Rbar", side = 2, line = 3)
lines(rbar_to_plot$mid.year, rbar_to_plot$rbar.eff, type = "l", col = COL_RBAR, lwd = LWD_RBAR)
abline(h = SSS_THRESHOLD, lty = "dotted", col = COL_SSS, lwd = 1.5)

# Añadir Sample Depth en el eje derecho
par(new = TRUE)
plot(txt_to_plot$year, txt_to_plot$samp.depth, type = "l", col = COL_SAMP_DEPTH, lwd = LWD_SAMP_DEPTH,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n",
     xlim = x_range_plot, ylim = c(0, max(txt_to_plot$samp.depth, na.rm=TRUE) * 1.1))
axis(4, lwd = 0, lwd.ticks = 1); mtext("Sample Depth", side = 4, line = 3, col = COL_SAMP_DEPTH)

# Eje X común
first_round_year_plot <- ceiling(START_PLOT_YEAR / 50) * 50
axis(1, at = seq(first_round_year_plot, max(txt_to_plot$year), by = 50), lwd = 0, lwd.ticks = 1)
mtext("Year", side = 1, line = 2.5)

### MODIFICADO: Estilo de la línea vertical ###
if (!is.na(first_year_threshold)) { abline(v = first_year_threshold, col = COL_SSS, lty = "dotted", lwd = 1.5) }
legend("topleft", 
       legend = c("SSS", paste0("Running Rbar (", WINDOW_RBAR, "-yr)"), "Sample Depth"), 
       col = c(COL_SSS, COL_RBAR, COL_SAMP_DEPTH), 
       lty = 1, lwd = 2, bty = "n", cex = 0.9)

# Rectángulo exterior
par(fig = c(0, 1, 0, 1), new = TRUE)
plot(0, 0, type = "n", xlab = "", ylab = "", axes = FALSE)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], border = "black", lwd = 1)

dev.off()
print(paste("Gráfico final con proporciones ajustadas guardado en:", file.path(OUTPUT_PATH, OUTPUT_FILE)))
