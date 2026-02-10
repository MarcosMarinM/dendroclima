# ==============================================================================
# 02_extract_pithoffsets.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Pith Offset Extractor for Age-Dependent Analysis.
#   Parses dendrochronological data files (specifically Besançon/TRDAS format) 
#   to extract "pith offset" estimates (the estimated number of rings missing 
#   to the pith).
#
#   Why this is important:
#   - Essential for Regional Curve Standardisation (RCS) and other age-dependent 
#     detrending methods.
#   - Allows calculation of the true biological age of each tree ring.
#
#   Inputs:
#   - Text file in Besançon format containing metadata tags (KeyCode, YearsToPith).
#
#   Outputs:
#   - CSV file mapping Series ID to Pith Offset values.
# ==============================================================================

# --- PARÁMETROS A MODIFICAR ---
DIR_TRABAJO      <- "PLACEHOLDER/path/to/working_dir"
INPUT_FILE_NAME  <- "PLACEHOLDER/path/to/input.txt"
OUTPUT_FILE_NAME <- "PLACEHOLDER/path/to/output.txt"
# ------------------------------


# Establecer el directorio de trabajo
setwd(DIR_TRABAJO)

# Función para leer el archivo y extraer la información necesaria
parse_file <- function(file_path) {
  lines <- readLines(file_path)
  series_list <- c()
  pithoffset_list <- c()
  pithoffset <- 1  # Valor por defecto
  
  for (line in lines) {
    if (grepl("KeyCode=", line)) {
      series <- sub(".*KeyCode=", "", line)
      series_list <- c(series_list, series)
      pithoffset_list <- c(pithoffset_list, pithoffset)  # Añadir el valor anterior de pithoffset
      pithoffset <- 1  # Reiniciar pithoffset al valor por defecto
    }
    if (grepl("YearsToPith=", line)) {
      pithoffset <- as.numeric(sub(".*YearsToPith=", "", line))
    }
  }
  
  # Eliminar el primer valor de pithoffset_list que es redundante
  pithoffset_list <- pithoffset_list[-1]
  
  # Asegurarse de que ambos vectores tengan la misma longitud
  if (length(series_list) != length(pithoffset_list)) {
    pithoffset_list <- c(pithoffset_list, pithoffset)
  }
  
  # Sumar 1 a todos los pithoffset excepto a los que son 1 porque están vacíos
  pithoffset_list <- sapply(pithoffset_list, function(po) if (po != 1) po + 1 else po)
  
  return(data.frame(series = series_list, pithoffset = pithoffset_list))
}

# Leer el archivo y obtener los datos
data <- parse_file(INPUT_FILE_NAME)

# Escribir los datos en un archivo CSV
write.csv(data, OUTPUT_FILE_NAME, row.names = FALSE)

cat("El archivo CSV ha sido creado en", file.path(DIR_TRABAJO, OUTPUT_FILE_NAME), "\n")
