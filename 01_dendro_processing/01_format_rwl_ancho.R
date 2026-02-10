# ==============================================================================
# 01_format_rwl_ancho.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Climate Data Reshaping Tool. 
#   Converts linear climate data (daily or monthly with 'yyyymmdd' format) into 
#   a wide matrix format (Years x Months). This wide format is often required 
#   by specialised dendroclimatic software or for easy visualisation of 
#   monthly trends across years.
#
#   Methodology:
#   - Reads input text file using data.table for efficiency.
#   - Parses dates to extract Year and Month.
#   - Pivots data from long (Date, Value) to wide (Year, Jan, Feb... Dec).
#   - Rounds values for cleaner output.
#
#   Inputs:
#   - Text file with at least two columns: Date and Value.
#
#   Outputs:
#   - Tab-delimited text file with Year index and 12 monthly columns.
# ==============================================================================

# --- PARÁMETROS A MODIFICAR ---
INPUT_FILE  <- "../data/GENERAL/precips/eobs01_general.txt"
OUTPUT_FILE <- "../data/GENERAL/precips/precipsok/eobs01.txt"
# ------------------------------

# Carga de librerías
if (!require("data.table")) install.packages("data.table")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")

library(data.table)
library(dplyr)
library(tidyr)

# Función principal
convertir_formato <- function(input_path, output_path) {
  cat("Procesando: ", input_path, "\n")
  
  # Cargar los datos
  data <- fread(input_path, header = FALSE, col.names = c("date", "precipitation"))
  
  # Extraer el año y el mes de la fecha
  data <- data %>%
    mutate(year = substr(as.character(date), 1, 4),
           month = as.numeric(substr(as.character(date), 5, 6))) %>%
    select(-date)
  
  # Pivotar a formato ancho
  data_wide <- data %>%
    pivot_wider(names_from = month, values_from = precipitation, values_fill = list(precipitation = NA)) %>%
    arrange(year)
  
  # Redondear
  data_wide <- data_wide %>%
    mutate(across(-year, round, 1))
  
  # Guardar
  fwrite(data_wide, file = output_path, sep = "\t", row.names = FALSE, col.names = FALSE)
  cat("Guardado en: ", output_path, "\n")
}

# Ejecución
if (file.exists(INPUT_FILE)) {
  convertir_formato(INPUT_FILE, OUTPUT_FILE)
} else {
  warning(paste("El archivo de entrada no existe:", INPUT_FILE))
}
