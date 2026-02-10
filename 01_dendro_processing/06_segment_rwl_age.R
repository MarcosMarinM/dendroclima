# ==============================================================================
# 06_segment_rwl_age.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   RWL Segmentation by Tree Age/Generation.
#   Splits a single RWL file into multiple sub-files based on the start year (age) 
#   of the series. Useful for analyzing different generations of trees separately 
#   (e.g., Young vs. Old) or for specific age-dependent signal analysis.
#
#   Classification Logic (Configurable):
#   - Young (Group 1): Series starting after 1925.
#   - Middle (Group 2): Series starting between 1825 and 1924.
#   - Old (Group 3): Series starting before 1825.
#
#   Inputs:
#   - Folder containing RWL files.
#
#   Outputs:
#   - Segmented RWL files (_yng.rwl, _mid.rwl, _old.rwl) in Output Folder.
# ==============================================================================

# --- PARÁMETROS A MODIFICAR ---
CARPETA_ENTRADA_RWL <- "PLACEHOLDER/path/to/input_dir"
CARPETA_SALIDA      <- "PLACEHOLDER/path/to/output_dir"

# Sufijos
SUFIJO_GRUPO1 <- "_yng" # Joven
SUFIJO_GRUPO2 <- "_mid" # Medio
SUFIJO_GRUPO3 <- "_old" # Antiguo

# Umbrales
UMBRAL_GRUPO1_INICIO <- 1925 # Grupo 1: >= este año
UMBRAL_GRUPO2_INICIO <- 1825 # Grupo 2: >= este año
UMBRAL_GRUPO2_FIN    <- 1924 # Grupo 2: <= este año
# ------------------------------

# Cargar paquetes necesarios
library(dplR)
library(tools)




# --- INICIO DEL SCRIPT ---

# Validar la carpeta de entrada
if (!dir.exists(CARPETA_ENTRADA_RWL)) {
  stop(paste("Error: La carpeta de entrada especificada (", CARPETA_ENTRADA_RWL, ") no existe.",
             "Por favor, verifica la ruta en la sección de parámetros."))
}

# Crear la carpeta de salida si no existe
if (CARPETA_SALIDA != "" && CARPETA_SALIDA != "." && !dir.exists(CARPETA_SALIDA)) {
  cat(paste("La carpeta de salida especificada (", CARPETA_SALIDA, ") no existe. Intentando crearla...\n"))
  dir.create(CARPETA_SALIDA, recursive = TRUE, showWarnings = TRUE)
  if (!dir.exists(CARPETA_SALIDA)) {
    stop(paste("No se pudo crear la carpeta de salida:", CARPETA_SALIDA,
               "Por favor, créala manualmente o verifica los permisos."))
  } else {
    cat(paste("Carpeta de salida", CARPETA_SALIDA, "creada/verificada.\n"))
  }
}

# Obtener la lista de archivos .rwl en la carpeta de entrada
lista_archivos_rwl <- list.files(path = CARPETA_ENTRADA_RWL, pattern = "\\.[rR][wW][lL]$", full.names = TRUE)

if (length(lista_archivos_rwl) == 0) {
  stop(paste("No se encontraron archivos .rwl en la carpeta:", CARPETA_ENTRADA_RWL))
}

cat(paste("\nSe encontraron", length(lista_archivos_rwl), "archivos .rwl para procesar:\n"))
for(archivo_path in lista_archivos_rwl) {
  cat(paste("-", basename(archivo_path), "\n"))
}
cat("\n--- Iniciando procesamiento en bucle ---\n")

# Bucle para procesar cada archivo .rwl encontrado
for (archivo_rwl_actual_path in lista_archivos_rwl) {
  
  nombre_archivo_actual <- basename(archivo_rwl_actual_path)
  cat(paste("\n-------------------------------------------------\n"))
  cat(paste("Procesando archivo:", nombre_archivo_actual, "\n"))
  cat(paste("-------------------------------------------------\n"))
  
  # Leer el archivo .rwl actual
  datos_rwl_original <- tryCatch({
    read.rwl(archivo_rwl_actual_path)
  }, error = function(e) {
    cat(paste("ERROR al leer el archivo", nombre_archivo_actual, ":", e$message, "\nSaltando este archivo.\n"))
    return(NULL) # Devuelve NULL si hay un error
  })
  
  # Si hubo un error al leer, saltar al siguiente archivo
  if (is.null(datos_rwl_original)) {
    next
  }
  
  cat(paste("Se han leído", ncol(datos_rwl_original), "series cronológicas de", nombre_archivo_actual, ".\n"))
  
  stats_series <- rwl.stats(datos_rwl_original)
  primer_anio_por_serie <- stats_series$first
  nombres_de_series <- colnames(datos_rwl_original) 
  
  if (length(nombres_de_series) != length(primer_anio_por_serie)) {
    cat(paste("ADVERTENCIA para", nombre_archivo_actual, ": El número de nombres de columna (",
              length(nombres_de_series),
              ") no coincide con la longitud de 'primer_anio_por_serie' (",
              length(primer_anio_por_serie), ").",
              "Esto podría causar problemas en la segmentación de este archivo. Saltando este archivo.\n"))
    next 
  }
  
  # Lógica de segmentación
  condicion_grupo1 <- primer_anio_por_serie >= UMBRAL_GRUPO1_INICIO & !is.na(primer_anio_por_serie)
  series_grupo1_ids <- nombres_de_series[which(condicion_grupo1)]
  
  condicion_grupo2 <- primer_anio_por_serie >= UMBRAL_GRUPO2_INICIO &
    primer_anio_por_serie <= UMBRAL_GRUPO2_FIN &
    !is.na(primer_anio_por_serie)
  series_grupo2_ids <- nombres_de_series[which(condicion_grupo2)]
  
  condicion_grupo3 <- primer_anio_por_serie < UMBRAL_GRUPO2_INICIO & !is.na(primer_anio_por_serie)
  series_grupo3_ids <- nombres_de_series[which(condicion_grupo3)]
  
  # Crear los data frames para cada nuevo archivo .rwl
  if (length(series_grupo1_ids) > 0) {
    if (!all(series_grupo1_ids %in% colnames(datos_rwl_original))) {
      cat(paste("ERROR CRÍTICO para", nombre_archivo_actual, ": Algunos IDs del Grupo 1 no son nombres de columna válidos. Saltando escritura para este grupo.\n"))
      datos_rwl_grupo1 <- data.frame(row.names = rownames(datos_rwl_original)) 
    } else {
      datos_rwl_grupo1 <- datos_rwl_original[, series_grupo1_ids, drop = FALSE]
    }
  } else {
    datos_rwl_grupo1 <- data.frame(row.names = rownames(datos_rwl_original))
  }
  
  if (length(series_grupo2_ids) > 0) {
    if (!all(series_grupo2_ids %in% colnames(datos_rwl_original))) {
      cat(paste("ERROR CRÍTICO para", nombre_archivo_actual, ": Algunos IDs del Grupo 2 no son nombres de columna válidos. Saltando escritura para este grupo.\n"))
      datos_rwl_grupo2 <- data.frame(row.names = rownames(datos_rwl_original))
    } else {
      datos_rwl_grupo2 <- datos_rwl_original[, series_grupo2_ids, drop = FALSE]
    }
  } else {
    datos_rwl_grupo2 <- data.frame(row.names = rownames(datos_rwl_original))
  }
  
  if (length(series_grupo3_ids) > 0) {
    if (!all(series_grupo3_ids %in% colnames(datos_rwl_original))) {
      cat(paste("ERROR CRÍTICO para", nombre_archivo_actual, ": Algunos IDs del Grupo 3 no son nombres de columna válidos. Saltando escritura para este grupo.\n"))
      datos_rwl_grupo3 <- data.frame(row.names = rownames(datos_rwl_original))
    } else {
      datos_rwl_grupo3 <- datos_rwl_original[, series_grupo3_ids, drop = FALSE]
    }
  } else {
    datos_rwl_grupo3 <- data.frame(row.names = rownames(datos_rwl_original))
  }
  
  # Definir nombres para los nuevos archivos .rwl
  nombre_base_actual_sin_ext <- file_path_sans_ext(nombre_archivo_actual)
  
  archivo_salida_grupo1_nombre <- paste0(nombre_base_actual_sin_ext, SUFIJO_GRUPO1, ".rwl")
  archivo_salida_grupo2_nombre <- paste0(nombre_base_actual_sin_ext, SUFIJO_GRUPO2, ".rwl")
  archivo_salida_grupo3_nombre <- paste0(nombre_base_actual_sin_ext, SUFIJO_GRUPO3, ".rwl")
  
  archivo_salida_grupo1_completo <- file.path(CARPETA_SALIDA, archivo_salida_grupo1_nombre)
  archivo_salida_grupo2_completo <- file.path(CARPETA_SALIDA, archivo_salida_grupo2_nombre)
  archivo_salida_grupo3_completo <- file.path(CARPETA_SALIDA, archivo_salida_grupo3_nombre)
  
  # Escribir los nuevos archivos .rwl
  if (ncol(datos_rwl_grupo1) > 0) {
    write.rwl(datos_rwl_grupo1, fname = archivo_salida_grupo1_completo, format = "tucson")
    cat(paste("Archivo creado:", archivo_salida_grupo1_completo, "con", ncol(datos_rwl_grupo1), "series.\n"))
  } else {
    cat(paste("Para", nombre_archivo_actual, ": No se encontraron series para el Grupo 1. No se creó el archivo.\n"))
  }
  
  if (ncol(datos_rwl_grupo2) > 0) {
    write.rwl(datos_rwl_grupo2, fname = archivo_salida_grupo2_completo, format = "tucson")
    cat(paste("Archivo creado:", archivo_salida_grupo2_completo, "con", ncol(datos_rwl_grupo2), "series.\n"))
  } else {
    cat(paste("Para", nombre_archivo_actual, ": No se encontraron series para el Grupo 2. No se creó el archivo.\n"))
  }
  
  if (ncol(datos_rwl_grupo3) > 0) {
    write.rwl(datos_rwl_grupo3, fname = archivo_salida_grupo3_completo, format = "tucson")
    cat(paste("Archivo creado:", archivo_salida_grupo3_completo, "con", ncol(datos_rwl_grupo3), "series.\n"))
  } else {
    cat(paste("Para", nombre_archivo_actual, ": No se encontraron series para el Grupo 3. No se creó el archivo.\n"))
  }
  
  # Resumen para el archivo actual
  cat(paste("\n--- Resumen para:", nombre_archivo_actual, "---\n"))
  cat(paste("Total de series en el archivo original:", ncol(datos_rwl_original), "\n"))
  cat(paste("Series en Grupo 1 (>= ", UMBRAL_GRUPO1_INICIO, "): ", length(series_grupo1_ids), "\n", sep=""))
  cat(paste("Series en Grupo 2 (", UMBRAL_GRUPO2_INICIO, "-", UMBRAL_GRUPO2_FIN, "): ", length(series_grupo2_ids), "\n", sep=""))
  cat(paste("Series en Grupo 3 (< ", UMBRAL_GRUPO2_INICIO, "): ", length(series_grupo3_ids), "\n", sep=""))
  
  total_series_asignadas <- length(series_grupo1_ids) + length(series_grupo2_ids) + length(series_grupo3_ids)
  series_con_na_en_primer_anio <- sum(is.na(primer_anio_por_serie))
  
  if (total_series_asignadas == (ncol(datos_rwl_original) - series_con_na_en_primer_anio)) {
    if (series_con_na_en_primer_anio > 0) {
      cat(paste(series_con_na_en_primer_anio, "serie(s) de", nombre_archivo_actual, "tenían un primer año NA y no fueron asignadas.\n"))
    }
  } else {
    cat(paste("ADVERTENCIA para", nombre_archivo_actual, ": El número de series asignadas (", total_series_asignadas,
              ") no coincide con el total original menos las series con primer año NA (",
              ncol(datos_rwl_original) - series_con_na_en_primer_anio,
              ").\n"))
  }
  cat(paste("--- Fin del procesamiento para:", nombre_archivo_actual, "---\n"))
  
} # Fin del bucle for

cat("\n\nProceso completado para todos los archivos.\n")