# ==============================================================================
# 14_reconstruct_nao.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Station-Based NAO Reconstructor.
#   A utility to build a specific NAO index reconstruction from two long instrumental 
#   pressure series (e.g., Azores and Reykjavik).
#
#   Process:
#   - Loads two raw pressure series.
#   - Standardizes each series relative to its own climatology (Monthly Z-scores).
#   - Calculates the difference (South Pole - North Pole).
#   - Re-standardizes the resulting difference to produce the final Index.
#
#   Inputs:
#   - Pressure file A (South Node, e.g., Azores).
#   - Pressure file B (North Node, e.g., Reykjavik).
#
#   Outputs:
#   - Reconstructed NAO Index (Year x Month matrix).
# ==============================================================================

library(dplyr)
library(tidyr)

# --- PARÁMETROS A MODIFICAR ---
INPUT_AZORES    <- "PLACEHOLDER/path/to/azores_pressure.txt"
INPUT_REIKIAVIK <- "PLACEHOLDER/path/to/reykjavik_pressure.txt"
OUTPUT_FILE     <- "PLACEHOLDER/path/to/nao_reconstructed.txt"
# ------------------------------


# 1. Leer datos (formato: 13 columnas sin cabecera)
azores <- read.table(INPUT_AZORES, header = FALSE)
reikiavik <- read.table(INPUT_REIKIAVIK, header = FALSE)

# 2. Función para normalizar por mes
procesar_datos <- function(df, nombre) {
  df %>%
    pivot_longer(
      cols = V2:V13,
      names_to = "mes_codigo",
      values_to = "presion"
    ) %>%
    mutate(
      mes = as.numeric(sub("V", "", mes_codigo)) - 1,
      año = V1
    ) %>%
    select(año, mes, presion) %>%
    group_by(mes) %>%
    mutate(
      media_clim = mean(presion, na.rm = TRUE),
      sd_clim = sd(presion, na.rm = TRUE),
      z_score = (presion - media_clim) / sd_clim  # Normalización
    ) %>%
    ungroup() %>%
    select(año, mes, !!sym(nombre) := z_score)
}

# 3. Procesar ambas estaciones
azores_z <- procesar_datos(azores, "z_azores")
reikiavik_z <- procesar_datos(reikiavik, "z_reikiavik")

# 4. Calcular NAO como diferencia y estandarizar
nao_index <- inner_join(
  azores_z,
  reikiavik_z,
  by = c("año", "mes")
) %>%
  mutate(
    dif_z = z_azores - z_reikiavik,
    nao = scale(dif_z)  # Estandarización final
  ) %>%
  select(año, mes, nao)

# 5. Pivotar los datos al formato deseado
nao_pivot <- nao_index %>%
  pivot_wider(names_from = mes, values_from = nao) %>%
  arrange(año)

# 6. Redondear los valores a tres decimales
nao_pivot <- nao_pivot %>%
  mutate(across(-año, ~ round(.x, 3)))

# 7. Exportación sin encabezado
write.table(nao_pivot, OUTPUT_FILE, sep = "\t", row.names = FALSE, col.names = FALSE)