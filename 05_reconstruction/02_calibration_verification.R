# ==============================================================================
# 02_calibration_verification.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Split-Period Calibration & Verification Analysis.
#   Validates the stability and predictive skill of the reconstruction model by 
#   splitting the instrumental period into two halves (Early / Late).
#
#   Procedure:
#   - Split-Sample Test:
#     1. Calibrate on Early Period -> Verify on Late Period.
#     2. Calibrate on Late Period -> Verify on Early Period.
#   - Metrics Calculated:
#     - R-squared (R2) for calibration and verification.
#     - Reduction of Error (RE).
#     - Coefficient of Efficiency (CE).
#
#   Visualisation:
#   - Generates a 4-panel plot showing Observed vs. Predicted values for both 
#     calibration and verification phases.
#
#   Inputs:
#   - Reconstructed Data File (containing Observed and Predicted columns).
#
#   Outputs:
#   - Calibration/Verification Plot (.pdf).
# ==============================================================================

# --- PARÁMETROS A MODIFICAR ---
INPUT_FILE_RECON <- 'PLACEHOLDER/path/to/reconstruction.txt'
OUTPUT_PDF_FILE  <- "PLACEHOLDER/path/to/plot.pdf"
# ------------------------------

# --- 2. Cargar Datos ---
# (Sin cambios, igual que antes)
tryCatch({
# --- 2. Cargar Datos ---
# (Sin cambios, igual que antes)
tryCatch({
  precip_data <- read.table(INPUT_FILE_RECON, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  precip_data$year <- as.numeric(precip_data$year)
  precip_data$res <- as.numeric(precip_data$res)
  precip_data$cumulative_precip <- as.numeric(precip_data$cumulative_precip)
  precip_data$reconstructed_precip <- as.numeric(precip_data$reconstructed_precip)
  print("Datos cargados correctamente.")
}, error = function(e) {
  stop("Error al cargar el archivo. Verifica la ruta y el separador. Error: ", e$message)
})

# ... (código intermedio sin cambios de variables de entrada) ...

# Abrir dispositivo PDF
pdf(OUTPUT_PDF_FILE, width = 8, height = 5, pointsize = base_pointsize, useDingbats = FALSE)

# Configurar ventana gráfica 2x2 y parámetros 
# mgp = c(label_pos, tick_num_pos, axis_line_pos)
par(mfrow = c(2, 2),
    mar = c(3.5, 3.5, 1.5, 1), # Margen interior reducido (abajo, izq, arriba, der)
    oma = c(0, 0, 0, 0),       # Sin margen exterior
    mgp = c(1.8, 0.6, 0),      # <-- ETIQUETAS (Year/Precip) MÁS PEGADAS AL EJE (valor inicial reducido)
    cex.axis = cex_axis_val,
    cex.lab = cex_lab_val
)

# --- Títulos formateados ---
r2_cal1_txt <- if (is.na(r2_cal1)) "NA" else sprintf("%.2f", r2_cal1)
r2_ver1_txt <- if (is.na(r2_ver1)) "NA" else sprintf("%.2f", r2_ver1)
r2_ver2_txt <- if (is.na(r2_ver2)) "NA" else sprintf("%.2f", r2_ver2)
r2_cal2_txt <- if (is.na(r2_cal2)) "NA" else sprintf("%.2f", r2_cal2)

title_cal1 <- bquote("Calibration (" * .(min(period1_years)) * "–" * .(max(period1_years)) * ")" * atop(R^2 == .(r2_cal1_txt)))
title_ver1 <- bquote("Verification (" * .(min(period2_years)) * "–" * .(max(period2_years)) * ")" * atop(R^2 == .(r2_ver1_txt)))
title_ver2 <- bquote("Verification (" * .(min(period1_years)) * "–" * .(max(period1_years)) * ")" * atop(R^2 == .(r2_ver2_txt)))
title_cal2 <- bquote("Calibration (" * .(min(period2_years)) * "–" * .(max(period2_years)) * ")" * atop(R^2 == .(r2_cal2_txt)))

# --- Posición para el texto dentro de los gráficos ---
# Usaremos coordenadas relativas (npc) convertidas a coordenadas de usuario
x_text_pos_npc <- 0.05 # 5% desde la izquierda
y_text_pos_npc <- 0.95 # 95% desde abajo (cerca de la parte superior)

# Gráfico 1: Calibración Ronda 1
plot(results_cal1$year, results_cal1$cumulative_precip, type = 'l', col = 'black',
     ylim = y_range, xlab = "Year", ylab = "Precipitation", axes = TRUE)
lines(results_cal1$year, results_cal1$predicted, col = 'red', lty = 2)
# Añadir título DENTRO con text()
x_pos <- grconvertX(x_text_pos_npc, from = "npc", to = "user")
y_pos <- grconvertY(y_text_pos_npc, from = "npc", to = "user")
text(x = x_pos, y = y_pos, labels = title_cal1, adj = c(0, 1), cex = cex_title_val) # adj=c(0,1) -> alinea esquina sup-izq del texto
legend("bottomright", legend = c("Observed", "Estimated"), col = c("black", "red"),
       lty = c(1, 2), bty = "n", cex = cex_legend_val) # Leyenda sólo aquí

# Gráfico 2: Verificación Ronda 1
plot(results_ver1$year, results_ver1$cumulative_precip, type = 'l', col = 'black',
     ylim = y_range, xlab = "Year", ylab = "Precipitation", axes = TRUE)
lines(results_ver1$year, results_ver1$predicted, col = 'red', lty = 2)
x_pos <- grconvertX(x_text_pos_npc, from = "npc", to = "user")
y_pos <- grconvertY(y_text_pos_npc, from = "npc", to = "user")
text(x = x_pos, y = y_pos, labels = title_ver1, adj = c(0, 1), cex = cex_title_val)

# Gráfico 3: Verificación Ronda 2
plot(results_ver2$year, results_ver2$cumulative_precip, type = 'l', col = 'black',
     ylim = y_range, xlab = "Year", ylab = "Precipitation", axes = TRUE)
lines(results_ver2$year, results_ver2$predicted, col = 'red', lty = 2)
x_pos <- grconvertX(x_text_pos_npc, from = "npc", to = "user")
y_pos <- grconvertY(y_text_pos_npc, from = "npc", to = "user")
text(x = x_pos, y = y_pos, labels = title_ver2, adj = c(0, 1), cex = cex_title_val)

# Gráfico 4: Calibración Ronda 2
plot(results_cal2$year, results_cal2$cumulative_precip, type = 'l', col = 'black',
     ylim = y_range, xlab = "Year", ylab = "Precipitation", axes = TRUE)
lines(results_cal2$year, results_cal2$predicted, col = 'red', lty = 2)
x_pos <- grconvertX(x_text_pos_npc, from = "npc", to = "user")
y_pos <- grconvertY(y_text_pos_npc, from = "npc", to = "user")
text(x = x_pos, y = y_pos, labels = title_cal2, adj = c(0, 1), cex = cex_title_val)


# Cerrar dispositivo PDF
dev.off()

# Resetear parámetros gráficos
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0),
    mgp = c(3, 1, 0), cex.axis = 1, cex.lab = 1, cex.main = 1.2, cex = 1)

cat(paste0("\nProceso completado. Gráfico v3 guardado como '", OUTPUT_PDF_FILE, "'.\n"))
