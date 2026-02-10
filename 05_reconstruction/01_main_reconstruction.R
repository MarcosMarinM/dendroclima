# ==============================================================================
# 01_main_reconstruction.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Annual Precipitation Reconstruction (Oct-Sep) with SSS Integration.
#   Performs the core step of reconstructing a specific climate variable (Annual 
#   Precipitation, Hydrological Year Oct-Sep) from a tree-ring chronology.
#
#   Methodology:
#   1. Climate Formatting: Aggregates monthly climate data into a Hydrological Year 
#      (Oct previous year to Sep current year).
#   2. Calibration: Fits a Linear Regression Model (Precip ~ RWI_res) over a 
#      specified calibration period (e.g., 1931-2019).
#   3. Reconstruction: Applies the model to the full length of the tree-ring chronology.
#   4. Uncertainty & Quality: Calculates SSS (Subsample Signal Strength) to assess 
#      reliability throughout time.
#   5. Visualisation: plots the reconstruction with 11-year smoothing and error ribbons.
#
#   Inputs:
#   - Monthly Climate Data file.
#   - Residual Chronology file (.txt and .rwi).
#
#   Outputs:
#   - Reconstruction Data File (.txt).
#   - Summary Plot (.pdf).
#   - Model Statistics and Equation (.txt).
# ==============================================================================

library(dplR)
library(dplyr)
library(ggplot2)
library(zoo)


# --- PARÁMETROS A MODIFICAR ---
# --- PARÁMETROS A MODIFICAR ---
INPUT_CLIMA_MENSUAL <- "PLACEHOLDER/path/to/climate_data.txt"
INPUT_CRONO_RES     <- "PLACEHOLDER/path/to/chronology.txt"
INPUT_RWI_RES       <- "PLACEHOLDER/path/to/residuals.rwi"
OUTPUT_FILE_TXT     <- "PLACEHOLDER/path/to/reconstruction.txt"
# ------------------------------

if (!dir.exists(dirname(OUTPUT_FILE_TXT))) dir.create(dirname(OUTPUT_FILE_TXT), recursive = TRUE, showWarnings = FALSE)


# ----------------------------------------------------------------------
# 1. CLIMA MENSUAL → SERIE OCT–SEP
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# 1. CLIMA MENSUAL → SERIE OCT–SEP
# ----------------------------------------------------------------------
clima <- read.table(INPUT_CLIMA_MENSUAL, header = TRUE, sep = "\t")

vars_cols <- names(clima)
col_variable <- vars_cols[!vars_cols %in% c("date")][1]
clima$valor <- clima[[col_variable]]

clima$year  <- as.numeric(substr(as.character(clima$date), 1, 4))
clima$month <- as.numeric(substr(as.character(clima$date), 5, 6))

# Año hidrológico: OCT (año-1) – SEP (año)
clima <- clima %>%
  mutate(
    hydro_year = ifelse(month >= 10, year + 1, year)
  )

clima_agg <- clima %>%
  group_by(hydro_year) %>%
  summarise(
    P_OctSep = sum(valor, na.rm = TRUE),
    .groups = "drop"
  )

# Periodo de calibración: 1931–2019
clima_cal <- clima_agg %>%
  filter(hydro_year >= 1931, hydro_year <= 2019)

# ----------------------------------------------------------------------
# 2. CRONOLOGÍA RESIDUAL (columna res)
# ----------------------------------------------------------------------
crono <- read.table(INPUT_CRONO_RES, header = TRUE, sep = "\t")
if (!"Year" %in% names(crono)) names(crono)[1] <- "Year"
if (!"res" %in% names(crono)) stop("La cronología no tiene columna 'res'.")

crono_res <- crono %>%
  select(Year, res) %>%
  mutate(Year = as.numeric(Year))

# Enlazamos Year con hydro_year (res de año N → clima Oct N-1–Sep N)
datos_cal <- inner_join(
  clima_cal,
  crono_res,
  by = c("hydro_year" = "Year")
)

if (nrow(datos_cal) < 10) stop("Muy pocos años en común para calibración.")

# ----------------------------------------------------------------------
# 3. MODELO DE CALIBRACIÓN Y RECONSTRUCCIÓN COMPLETA
# ----------------------------------------------------------------------
mod <- lm(P_OctSep ~ res, data = datos_cal)

# Reconstrucción para todos los años donde haya res
crono_res_todo <- crono_res %>%
  rename(hydro_year = Year)

recon <- crono_res_todo %>%
  mutate(
    P_recon = predict(mod, newdata = data.frame(res = res))
  )

# ----------------------------------------------------------------------
# 4. SSS AÑO A AÑO DESDE RWI RESIDUAL
# ----------------------------------------------------------------------
rwi <- read.rwl(INPUT_RWI_RES)

sss_vec <- sss(rwi)  # vector con nombres = años
years_sss <- as.numeric(names(sss_vec))

sss_df <- data.frame(
  hydro_year = years_sss,
  SSS = as.numeric(sss_vec)
)

# ----------------------------------------------------------------------
# 5. UNIR RECONSTRUCCIÓN + SSS Y GUARDAR
# ----------------------------------------------------------------------
out <- recon %>%
  left_join(sss_df, by = "hydro_year") %>%
  arrange(hydro_year)

# Guardar tabla principal
write.table(
  out,
  file = OUTPUT_FILE_TXT,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ======================================================================
# 6. GRÁFICO EXACTO ESTILO TU CÓDIGO (solo SSS >= 0.85)
# ======================================================================


# Definir out_full para que out_plot funcione (parece que falta en el original, asumo out es out_full referenciado)
out_full <- out 

out_plot <- out_full %>%
  filter(!is.na(SSS), SSS >= 0.85)

# Calcular media móvil 11 años de reconstrucción
out_plot$media_movil_11 <- rollmean(out_plot$P_recon, k = 11, fill = NA, align = "center")

# RMSE sobre periodo calibración (donde hay observado)
rmse <- sqrt(mean((out_plot$P_OctSep - out_plot$P_recon)^2, na.rm = TRUE))

# Años extremos por siglo (más bajo y más alto)
out_plot$Century <- (out_plot$hydro_year %/% 100) * 100
lowest_years <- out_plot %>% 
  group_by(Century) %>% 
  filter(P_recon == min(P_recon, na.rm = TRUE)) %>%
  slice(1)  # solo el primero si hay empates

highest_years <- out_plot %>% 
  group_by(Century) %>% 
  filter(P_recon == max(P_recon, na.rm = TRUE)) %>%
  slice(1)

# Gráfico EXACTO como tu código
p <- ggplot(out_plot, aes(x = hydro_year)) +
  geom_line(aes(y = P_recon, color = "Predicted"), linewidth = 0.4) +
  geom_line(aes(y = P_OctSep, color = "Observed"), linewidth = 0.4, linetype = "dashed") +
  geom_line(aes(y = media_movil_11), color = "#8B6969", linewidth = 0.5) +
  geom_hline(yintercept = mean(out_plot$P_OctSep, na.rm = TRUE), linetype = "dashed", color = "#CDAD00") +
  geom_ribbon(aes(ymin = P_recon - rmse, ymax = P_recon + rmse), fill = "#999", alpha = 0.4) +
  
  # Etiquetas años extremos por siglo
  geom_text(data = lowest_years, aes(x = hydro_year, y = P_recon, label = hydro_year), 
            angle = 45, vjust = 1.5, hjust = 1, size = 3, colour = "red") +
  geom_text(data = highest_years, aes(x = hydro_year, y = P_recon, label = hydro_year), 
            angle = 45, vjust = -0.5, hjust = 0, size = 3, colour = "blue") +
  
  scale_color_manual(values = c("Observado" = "red", "Predicho" = "blue")) +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  labs(
    x = "Hydrological year (Oct-Sep)",
    y = "Precipitation (mm)",
    color = "Legend",
    title = "Precipitation reconstruction since 931 CE"
  ) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid"))))

ruta_plot <- gsub("\\.txt$", "_Grafico_SSS085.pdf", OUTPUT_FILE_TXT)
ggsave(p, filename = ruta_plot, width = 12, height = 6)

# ======================================================================
# 7. GUARDAR ECUACIÓN DEL MODELO (igual que antes)
# ======================================================================
coef_mod <- coef(mod)
eq_str <- sprintf("P_OctSep = %.4f + %.4f * res", coef_mod[1], coef_mod[2])

ruta_eq <- gsub("\\.txt$", "_EcuacionModelo.txt", OUTPUT_FILE_TXT)
cat(eq_str, "\nR² =", round(summary(mod)$r.squared, 4), "\nRMSE =", round(rmse, 2), "\n\n", 
    file = ruta_eq)
capture.output(summary(mod), file = ruta_eq, append = TRUE)

cat("\nArchivos generados:\n")
cat("- Tabla:", OUTPUT_FILE_TXT, "\n")
cat("- Gráfico:", ruta_plot, "\n")
cat("- Ecuación:", ruta_eq, "\n")
