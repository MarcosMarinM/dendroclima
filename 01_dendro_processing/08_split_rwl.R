# ==============================================================================
# 08_split_rwl.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Macro-Collection Splitter.
#   Deconstructs a massive "Macro-RWL" file (containing series from many different sites) 
#   into individual RWL files per site.
#
#   Mechanism:
#   - Uses a metadata table to map Site IDs (prefixes) to their corresponding series.
#   - Filter columns in the master RWL based on these prefixes.
#   - Exports clean, site-specific RWL files.
#
#   Inputs:
#   - Macro-RWL file (all series).
#   - Metadata text file (Collection definition).
#
#   Outputs:
#   - Individual .rwl files for each site in the collection.
# ==============================================================================

# --- PARÁMETROS A MODIFICAR ---
PATH_MACRO_RWL <- "PLACEHOLDER/path/to/macro_rwl.rwl"
PATH_METADATA  <- "PLACEHOLDER/path/to/metadata_collection.txt"
DIR_OUTPUT     <- "PLACEHOLDER/path/to/output_dir"
# ------------------------------


# ------------------------------------------------------------------------------
# 2. CARGA DE LIBRERÍAS
# ------------------------------------------------------------------------------
library(dplR)
library(dplyr)
library(stringr)

# ------------------------------------------------------------------------------
# 3. LECTURA DE DATOS
# ------------------------------------------------------------------------------

cat(">>> Cargando datos...\n")

# Cargar la macrocolección
if (file.exists(PATH_MACRO_RWL)) {
  macro_rwl <- read.rwl(PATH_MACRO_RWL)
  cat("    Macrocolección cargada: ", ncol(macro_rwl), " series encontradas.\n")
} else {
  stop("ERROR: No se encuentra el archivo RWL en la ruta especificada.")
}

# Cargar los metadatos
if (file.exists(PATH_METADATA)) {
  meta <- read.table(PATH_METADATA, header = TRUE, stringsAsFactors = FALSE)
  cat("    Metadatos cargados: ", nrow(meta), " sitios definidos.\n")
} else {
  stop("ERROR: No se encuentra el archivo de metadatos.")
}

# Crear directorio de salida
if (!dir.exists(DIR_OUTPUT)) {
  dir.create(DIR_OUTPUT, recursive = TRUE)
  cat("    Directorio de salida creado: ", DIR_OUTPUT, "\n")
}

cat("------------------------------------------------------------------------------\n")

# ------------------------------------------------------------------------------
# 4. PROCESAMIENTO Y GENERACIÓN DE ARCHIVOS
# ------------------------------------------------------------------------------

# Obtenemos la lista única de colecciones originales (columna 'rwl' y su prefijo 'pre')
lista_sitios <- meta %>%
  select(rwl, pre) %>%
  distinct()

cat(">>> Iniciando la separación de series...\n\n")

# Bucle para procesar cada sitio
for (i in 1:nrow(lista_sitios)) {
  
  nombre_archivo_destino <- lista_sitios$rwl[i]
  prefijo_busqueda <- lista_sitios$pre[i]
  
  # 4.1. Identificar las series
  patron <- paste0("^", prefijo_busqueda)
  nombres_series_rwl <- colnames(macro_rwl)
  indices_coincidencia <- grep(pattern = patron, x = nombres_series_rwl, ignore.case = TRUE)
  
  # 4.2. Verificar si se encontraron series
  if (length(indices_coincidencia) > 0) {
    
    subset_rwl <- macro_rwl[, indices_coincidencia]
    subset_rwl <- subset_rwl[rowSums(!is.na(subset_rwl)) > 0, ]
    
    # 4.3. Escribir el nuevo archivo .rwl
    nombre_final <- if(str_ends(nombre_archivo_destino, "\\.rwl")) nombre_archivo_destino else paste0(nombre_archivo_destino, ".rwl")
    ruta_guardado <- file.path(DIR_OUTPUT, nombre_final)
    
    write.rwl(rwl.df = subset_rwl, fname = ruta_guardado, format = "tucson")

    
    # Log de progreso
    cat(sprintf("[OK] %s guardado | Prefijo: '%s' | Series: %d\n", 
                nombre_final, prefijo_busqueda, ncol(subset_rwl)))
    
    archivos_creados <- archivos_creados + 1
    series_totales_exportadas <- series_totales_exportadas + ncol(subset_rwl)
    
  } else {
    # Alerta si hay un sitio en la lista pero no series en el RWL
    cat(sprintf("[ALERTA] No se encontraron series para '%s' (Prefijo: '%s')\n", 
                nombre_archivo_destino, prefijo_busqueda))
  }
}

# ------------------------------------------------------------------------------
# 5. RESUMEN FINAL
# ------------------------------------------------------------------------------
cat("\n------------------------------------------------------------------------------\n")
cat("PROCESO TERMINADO.\n")
cat("Resumen:\n")
cat(" - Archivos creados:", archivos_creados, "\n")
cat(" - Total de series procesadas:", series_totales_exportadas, " de ", ncol(macro_rwl), "\n")
cat(" - Ubicación:", DIR_OUTPUT, "\n")


# Verificación de integridad simple
if (series_totales_exportadas < ncol(macro_rwl)) {
  cat("\nNOTA: El número de series exportadas es menor al total original.\n")
  cat("Esto puede ser normal si hay series en la macrocolección sin metadatos asociados\n")
  cat("o si los prefijos no coinciden exactamente.\n")
}