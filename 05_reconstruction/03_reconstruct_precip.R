# ==============================================================================
# 03_reconstruct_precip.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Precision Daily-Window Precipitation Reconstruction.
#   An advanced reconstruction script that calculates cumulative precipitation 
#   from DAILY data for a specific, flexible seasonal window (e.g., 320 days 
#   ending on June 30th) before building the reconstruction model.
#
#   Features:
#   1. Flexible Window: Computes cumulative precip for any custom window defined 
#      by end date and length in days.
#   2. SSS Calculation: Calculates Subsample Signal Strength directly from the .rwl file.
#   3. Full Pipeline:
#      - Aggregates daily precip -> Annual values.
#      - Merges with Tree-Ring Data.
#      - Calibrates Linear Model.
#      - Reconstructs full series.
#      - Combines with SSS.
#      - Visualizes results (Observed vs Reconstructed + Moving Average).
#
#   Inputs:
#   - Daily Precipitation File.
#   - Tree-Ring Data File (.txt).
#   - Raw RWL File (for SSS).
#
#   Outputs:
#   - Detailed Reconstruction File (.txt).
# ==============================================================================

# --- 0. Cargar/Instalar Librerías Necesarias ---

# --- 0. Cargar/Instalar Librerías Necesarias ---
library(dplyr)
library(lubridate)
library(ggplot2)
library(zoo)
library(dplR)

print("Librerías cargadas correctamente.")


# --- PARÁMETROS A MODIFICAR ---
# Ventana precipitación
WINDOW_DAYS      <- 320
WINDOW_END_MONTH <- 6  # Junio
WINDOW_END_DAY   <- 30 # Día 30

# Rutas
# Rutas
INPUT_PRECIP_FILE    <- 'PLACEHOLDER/path/to/precip_daily.txt'
INPUT_TREE_RING_FILE <- 'PLACEHOLDER/path/to/tree_ring_data.txt'
INPUT_CRONO_RWL      <- "PLACEHOLDER/path/to/chronology.rwl"
OUTPUT_FILE          <- 'PLACEHOLDER/path/to/reconstruction_output.txt'
# ------------------------------


# --- 2. Cargar datos de precipitación diaria ---

precip_data <- read.table(INPUT_PRECIP_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
precip_data$date <- ymd(precip_data$date) # Convertir la columna 'date' a formato Fecha

print("Primeras filas de datos de precipitación:")
head(precip_data)
print("Rango de fechas de precipitación:")
range(precip_data$date, na.rm = TRUE)

# --- 3. Calcular precipitación acumulada anual ---

print(paste("Calculando precipitación acumulada para una ventana de", WINDOW_DAYS,
            "días terminando el", WINDOW_END_DAY, "de", month.name[WINDOW_END_MONTH]))

annual_precip_sum <- data.frame(year = integer(),
                                cumulative_precip = numeric())

min_date_available <- min(precip_data$date, na.rm = TRUE)
max_date_available <- max(precip_data$date, na.rm = TRUE)

# Determinar los años posibles para el cálculo
# Se necesita empezar un año después del primer dato para tener un período completo previo
start_year_possible <- year(min_date_available) + 1
end_year_possible <- year(max_date_available)

# Ajustar el último año si no hay datos hasta la fecha final de la ventana
last_potential_end_date <- ymd(sprintf("%d-%02d-%02d", end_year_possible, WINDOW_END_MONTH, WINDOW_END_DAY))
if (max_date_available < last_potential_end_date) {
  end_year_possible <- end_year_possible - 1
}

print(paste("Años potenciales para cálculo de precipitación acumulada:", start_year_possible, "a", end_year_possible))

# Calcular fechas de inicio/fin genéricas (solo para el subtítulo del gráfico)
generic_end_date_for_subtitle <- ymd(sprintf("2000-%02d-%02d", WINDOW_END_MONTH, WINDOW_END_DAY)) # Año irrelevante
generic_start_date_for_subtitle <- generic_end_date_for_subtitle - days(WINDOW_DAYS - 1)


for (current_year in start_year_possible:end_year_possible) {
  # Calcular fecha de fin y de inicio para el año actual
  end_date <- ymd(sprintf("%d-%02d-%02d", current_year, WINDOW_END_MONTH, WINDOW_END_DAY))
  # La fecha de inicio es N-1 días antes de la fecha final para una ventana de N días
  start_date <- end_date - days(WINDOW_DAYS - 1)
  
  # Procesar solo si el período completo está dentro de los datos disponibles
  if (start_date >= min_date_available && end_date <= max_date_available) {
    period_data <- precip_data %>%
      filter(date >= start_date & date <= end_date)
    
    # Asegurarse de que realmente hay datos cubriendo los extremos del período
    if (nrow(period_data) > 0 && min(period_data$date, na.rm = TRUE) <= start_date && max(period_data$date, na.rm = TRUE) >= end_date) {
      total_precip <- sum(period_data$precipitation, na.rm = TRUE)
      annual_precip_sum <- rbind(annual_precip_sum,
                                 data.frame(year = current_year,
                                            cumulative_precip = total_precip))
    } else {
      # Advertir si faltan datos DENTRO del período (huecos, etc.)
      print(paste("Datos incompletos para el período completo del año:", current_year, "(Fechas:", start_date, "a", end_date, ")"))
    }
  } 
}

print("Precipitación acumulada anual calculada (primeras filas):")
head(annual_precip_sum)
print("Años con precipitación acumulada calculada:")
if (nrow(annual_precip_sum) > 0) print(sort(annual_precip_sum$year)) else print("Ninguno")


# --- 4. Cargar datos de anillos (txt) ---

tree_ring_data_full <- read.table(INPUT_TREE_RING_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Seleccionar solo las columnas necesarias y quitar NAs en 'res'
tree_ring_data <- tree_ring_data_full %>%
  select(year, std, res) %>%
  filter(!is.na(res))

print("Primeras filas de datos de anillos (std y res):")
head(tree_ring_data)
print("Rango de años en los datos de anillos:")
range(tree_ring_data$year, na.rm = TRUE)

# --- 5. Unir datos para calibración ---

# Combinar precipitación calculada y datos de anillos por año
calibration_data <- inner_join(annual_precip_sum, tree_ring_data, by = "year")

print("Datos combinados para calibración (primeras filas):")
head(calibration_data)
print(paste("Número de años disponibles para calibración:", nrow(calibration_data)))
if (nrow(calibration_data) > 0) print(paste("Rango de años de calibración:", min(calibration_data$year, na.rm = TRUE), "-", max(calibration_data$year, na.rm = TRUE)))


# --- 6. Modelo de regresión lineal (calibración) ---

# Necesitamos al menos 2 puntos para una regresión lineal
if (nrow(calibration_data) < 2) {
  stop("No hay suficientes años solapados (>= 2) para calibrar el modelo.")
}

# Modelo: Predecir precipitación acumulada a partir del índice residual 'res'
calibration_model <- lm(cumulative_precip ~ res, data = calibration_data)

print("Resumen del modelo de calibración:")
summary_model <- summary(calibration_model)
print(summary_model)

# Extraer métricas clave del modelo
r_squared <- summary_model$r.squared
adj_r_squared <- summary_model$adj.r.squared
if (nrow(summary_model$coefficients) > 1) {
  p_value_model <- summary_model$coefficients[2, 4] # P-valor del predictor 'res'
  print(paste("P-valor del predictor (res):", format.pval(p_value_model, digits = 3)))
} else {
  print("Advertencia: El modelo solo tiene intercepto, no hay predictor 'res'.")
  p_value_model <- NA
}
print(paste("R-cuadrado:", round(r_squared, 3)))
print(paste("R-cuadrado ajustado:", round(adj_r_squared, 3)))


# --- 7. Reconstruir precipitación ---

# Usar el modelo calibrado para predecir sobre todos los datos de anillos
reconstruction <- tree_ring_data %>%
  mutate(reconstructed_precip = predict(calibration_model, newdata = .))

# --- 8. Calcular SSS (Sample Size Statistic) ---
print(paste("Cargando archivo RWL para cálculo de SSS:", INPUT_CRONO_RWL))

sss_df <- data.frame(year = integer(), sss = numeric()) # Inicializar vacío
sss_calculated <- FALSE

if (file.exists(INPUT_CRONO_RWL)) {
  tryCatch({
    cronologia_rwl <- read.rwl(INPUT_CRONO_RWL)
    sss_values <- sss(cronologia_rwl) # Calcula el número de series por año
    
    # Convertir el resultado (vector nombrado) a data frame
    years_from_sss <- suppressWarnings(as.numeric(names(sss_values)))
    valid_years_indices <- !is.na(years_from_sss)
    
    if (length(sss_values) > 0 && any(valid_years_indices)) {
      sss_df <- data.frame(
        year = years_from_sss[valid_years_indices],
        sss = as.numeric(sss_values[valid_years_indices])
      )
      print("SSS calculado y convertido a data frame (primeras filas):")
      print(head(sss_df))
      print(paste("Rango de años con SSS:", min(sss_df$year, na.rm = T), "-", max(sss_df$year, na.rm = T)))
      sss_calculated <- TRUE
    } else {
      print("ADVERTENCIA: No se pudieron extraer años válidos o el vector SSS está vacío.")
    }
  }, error = function(e) {
    # Capturar error si falla read.rwl o sss
    print(paste("ERROR al leer RWL o calcular SSS:", e$message))
  })
} else {
  print(paste("¡ADVERTENCIA! No se encontró el archivo RWL:", INPUT_CRONO_RWL))
  print("La columna 'sss' quedará vacía (NA).")
}


# --- 9. Combinar resultados y calcular media móvil ---

# Unir datos de anillos, precipitación instrumental y reconstrucción
final_results <- full_join(tree_ring_data, annual_precip_sum, by = "year")
final_results <- full_join(final_results, select(reconstruction, year, reconstructed_precip), by = "year")

# Añadir la columna SSS (si se calculó)
if (sss_calculated && nrow(sss_df) > 0) {
  final_results <- left_join(final_results, sss_df, by = "year")
} else {
  # Añadir columna SSS con NAs si no se pudo calcular
  final_results <- final_results %>% mutate(sss = NA_real_)
}

# Ordenar por año y calcular media móvil centrada de 11 años
final_results <- final_results %>%
  arrange(year) %>%
  mutate(reconstructed_precip_mov_avg_11yr = rollmean(reconstructed_precip, k = 11, fill = NA, align = "center"))

print("Resultados finales combinados (primeras filas):")
head(final_results)
print("Resultados finales combinados (últimas filas):")
tail(final_results)
print("Nombres de las columnas finales:")
print(names(final_results))

# --- 10. Guardar resultados ---

# Seleccionar y ordenar columnas para el archivo de salida
output_data <- final_results %>%
  select(year, res, std, sss, cumulative_precip, reconstructed_precip, reconstructed_precip_mov_avg_11yr)

# Guardar como archivo de texto delimitado por tabuladores
write.table(output_data,
            file = OUTPUT_FILE,
            sep = "\t",
            row.names = FALSE, # No incluir números de fila
            quote = FALSE,     # No usar comillas
            na = "NA")         # Representar NAs como "NA"

print(paste("Resultados guardados en:", OUTPUT_FILE))

# --- 11. Visualización ---

# Gráfico principal: Reconstrucción vs Instrumental
plot_reconstruction <- ggplot(final_results, aes(x = year)) +
  geom_line(aes(y = reconstructed_precip, color = "Reconstrucción"), size = 0.8) +
  geom_line(aes(y = reconstructed_precip_mov_avg_11yr, color = "Media móvil 11 años"), size = 1.2, linetype = "dashed") +
  geom_point(aes(y = cumulative_precip, color = "Instrumental"), size = 0.5, alpha = 0.8) + # Puntos más pequeños
  scale_color_manual(
    name = "Serie",
    values = c("Reconstrucción" = "steelblue2", "Instrumental" = "firebrick", "Media móvil 11 años" = "grey50"),
    breaks = c("Instrumental", "Reconstrucción", "Media móvil 11 años") # Orden leyenda
  ) +
  labs(
    title = "Reconstrucción dendroclimatológica de precipitación acumulada",
    subtitle = paste("Período:", format(generic_start_date_for_subtitle, "%d %b"),
                     "a", format(generic_end_date_for_subtitle, "%d %b"),
                     paste0("(", WINDOW_DAYS, " días)")),
    x = "Año",
    y = "Precipitación acumulada (mm)",
    caption = if (nrow(calibration_data) > 0) {
      paste("Calibración:", min(calibration_data$year, na.rm = T), "-", max(calibration_data$year, na.rm = T),
            " R²adj:", round(adj_r_squared, 2), " p:", format.pval(p_value_model, digits = 2))
    } else {
      "Sin datos suficientes para calibración"
    }
  ) +
  theme_bw() + # Tema limpio
  theme(legend.position = "bottom") # Leyenda abajo

print(plot_reconstruction)

# Gráfico opcional SSS
if (sss_calculated && nrow(final_results %>% filter(!is.na(sss))) > 0) {
  plot_sss <- ggplot(final_results %>% filter(!is.na(sss)), aes(x = year, y = sss)) +
    geom_line(color = "darkgreen") +
    # geom_point(color = "darkgreen", size = 1.5) + # Opcional: añadir puntos
    labs(
      title = "SSS", # Título más descriptivo
      x = "Año",
      y = "SSS"
    ) +
    theme_bw()
  print(plot_sss)
} else {
  print("No se generó gráfico SSS (sin datos válidos).")
}

# Guardar gráficos (opcional)
# ggsave("reconstruccion_precipitacion.png", plot = plot_reconstruction, width = 11, height = 7)
# if (exists("plot_sss")) ggsave("sss_evolucion.png", plot = plot_sss, width = 10, height = 5)


print("--- FIN DEL SCRIPT ---")