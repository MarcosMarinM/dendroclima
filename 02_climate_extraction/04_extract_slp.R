# ==============================================================================
# 05_extract_slp.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Sea Level Pressure (SLP) Data Formatter.
#   Proceses raw SLP data (often from online portals like KNMI or NOAA) into a 
#   standardised matrix format (Year x Month) for dendroclimatic analysis.
#
#   Steps:
#   - Reads raw text data.
#   - Selects the grid point nearest to a target coordinate (if multiple are present).
#   - Converts linear time steps to Calendar Year/Month.
#   - Reshapes data to Wide Format (Jan...Dec columns).
#
#   Inputs:
#   - Raw SLP text file.
#   - Target Lat/Lon.
#
#   Outputs:
#   - Formatted SLP matrix file.
# ==============================================================================

# --- PARÁMETROS A MODIFICAR ---
library(dplyr)
library(tidyr)

# --- PARÁMETROS A MODIFICAR ---
INPUT_FILE  <- "PLACEHOLDER/path/to/slp.txt"
OUTPUT_FILE <- "PLACEHOLDER/path/to/output.txt"
TARGET_LAT  <- 13.51
TARGET_LON  <- 2.1
# ------------------------------




# Leer los datos
datos <- read.table(INPUT_FILE, header = TRUE, sep = "\t")

# Calcular distancia a las coordenadas objetivo
datos$distancia <- sqrt((datos$lat - TARGET_LAT)^2 + (datos$lon - TARGET_LON)^2)

# Seleccionar el punto más cercano para cada tiempo
seleccionados <- datos %>%
  group_by(time) %>%
  slice(which.min(distancia)) %>%
  ungroup()

# Convertir tiempo a año y mes
seleccionados <- seleccionados %>%
  mutate(
    año = 1658 + floor((time + 11) / 12),
    mes = ((time + 11) %% 12) + 1
  )

# Crear todas las combinaciones posibles de año-mes
años <- min(seleccionados$año):max(seleccionados$año)
combinaciones <- expand.grid(año = años, mes = 1:12)

# Combinar con los datos seleccionados
datos_completos <- left_join(combinaciones, 
                             seleccionados %>% select(año, mes, slp),
                             by = c("año", "mes"))

# Reorganizar datos en formato ancho
tabla_final <- datos_completos %>%
  pivot_wider(
    names_from = mes,
    values_from = slp
  ) %>%
  arrange(año)

# Renombrar columnas de meses
colnames(tabla_final)[2:13] <- month.abb

# Reordenar columnas para diciembre como última
tabla_final <- tabla_final %>%
  select(año, Jan, Feb, Mar, Apr, May, Jun, 
         Jul, Aug, Sep, Oct, Nov, Dec)

# Escribir archivo de salida
write.table(tabla_final, OUTPUT_FILE, sep = "\t", 
            row.names = FALSE, col.names = FALSE, na = "NA")
