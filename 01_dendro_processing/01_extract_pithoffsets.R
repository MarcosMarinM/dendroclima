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

# Función para leer el archivo y extraer la información necesaria.
# Parsea bloques entre KeyCode tags, buscando YearsToPith dentro de cada bloque.
# Funciona independientemente del orden relativo de KeyCode y YearsToPith.
parse_file <- function(file_path) {
  lines <- readLines(file_path)
  
  # Find line indices where KeyCode appears
  key_lines <- grep("KeyCode=", lines)
  if (length(key_lines) == 0) stop("No se encontraron líneas con 'KeyCode=' en el archivo.")
  
  series_list <- character(length(key_lines))
  pithoffset_list <- numeric(length(key_lines))
  
  for (i in seq_along(key_lines)) {
    # Extract series name
    series_list[i] <- sub(".*KeyCode=", "", lines[key_lines[i]])
    
    # Define the block of lines belonging to this series
    start <- key_lines[i]
    end <- if (i < length(key_lines)) key_lines[i + 1] - 1 else length(lines)
    block <- lines[start:end]
    
    # Search for YearsToPith within the block
    ytp_line <- grep("YearsToPith=", block, value = TRUE)
    if (length(ytp_line) > 0) {
      po <- as.numeric(sub(".*YearsToPith=", "", ytp_line[1]))
      # Sumar 1 al offset (convención: YearsToPith=0 significa que se llegó al pith)
      pithoffset_list[i] <- po + 1
    } else {
      pithoffset_list[i] <- 1  # Valor por defecto si no hay YearsToPith
    }
  }
  
  data.frame(series = series_list, pithoffset = pithoffset_list)
}

# Leer el archivo y obtener los datos
data <- parse_file(INPUT_FILE_NAME)

# Escribir los datos en un archivo CSV
write.csv(data, OUTPUT_FILE_NAME, row.names = FALSE)

cat("El archivo CSV ha sido creado en", file.path(DIR_TRABAJO, OUTPUT_FILE_NAME), "\n")
