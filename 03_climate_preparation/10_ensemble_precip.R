# ==============================================================================
# 15_ensemble_precip.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Ensemble Precipitation Generator.
#   Combines multiple independent precipitation datasets (e.g., Station data, 
#   Gridded Reanalysis, Satellite products) to create a robust "Ensemble" series.
#
#   Methodology:
#   - Reads multiple input files.
#   - Aligns them by Year and Month.
#   - Computes the 'Ensemble Mean' and 'Ensemble Median' for each time step.
#   - Uses the Median to reduce the influence of outliers/extreme errors in single datasets.
#
#   Inputs:
#   - List of precipitation text files to combine.
#
#   Outputs:
#   - Ensemble Mean File.
#   - Ensemble Median File.
# ==============================================================================

# --- PARÁMETROS A MODIFICAR ---
CARPETA_INPUTS <- "PLACEHOLDER/path/to/precip_inputs"
ARCHIVOS_USAR  <- c("dataset1.txt", "dataset2.txt", "dataset3.txt")
# ------------------------------

# ===========================================================================
# Script para generar los archivos de precipitación "ensemble mean" y "ensemble median"
# ===========================================================================

# Función para generar archivos ensemble mean y median mensuales
generar_ensemble_mensual <- function(carpeta_inputs, archivos_usar) {

  # Inicializar lista para almacenar datos
  lista_datos <- list()
  
  # Leer cada archivo en la lista y añadirlo a la lista de datos
  for (archivo in archivos_usar) {
    ruta_archivo <- file.path(carpeta_inputs, archivo)
    if (file.exists(ruta_archivo)) {
      datos <- read.table(ruta_archivo, header = FALSE)
      colnames(datos) <- c("Year", paste0("Month_", 1:12))  # Asignar nombres a las columnas
      lista_datos[[archivo]] <- datos
    } else {
      warning(paste("El archivo", archivo, "no existe en la carpeta:", carpeta_inputs))
    }
  }
  
  # Comprobar si hay al menos un archivo cargado
  if (length(lista_datos) == 0) {
    stop("No se cargaron datos. Verifica los nombres de los archivos y la carpeta de inputs.")
  }
  
  # Combinar todos los datos en una sola tabla (permite años no comunes)
  datos_combinados <- Reduce(function(x, y) merge(x, y, by = "Year", all = TRUE), lista_datos)
  
  # Reorganizar los datos para tener columnas de cada mes en listas separadas
  datos_mensuales <- lapply(1:12, function(mes) {
    # Extraer columnas del mes actual de todos los datasets
    columnas_mes <- grep(paste0("Month_", mes), colnames(datos_combinados))
    datos_mes <- datos_combinados[, c(1, columnas_mes), drop = FALSE]  # Incluye columna Year y las del mes
    colnames(datos_mes)[-1] <- paste0("Dataset_", 1:(ncol(datos_mes) - 1))
    datos_mes
  })
  
  # Calcular medias y medianas mensuales
  resultados <- lapply(1:12, function(mes) {
    datos_mes <- datos_mensuales[[mes]]
    
    # Aplicar cálculos por fila (año) para el mes actual
    resultado_mes <- apply(datos_mes[, -1], 1, function(fila) {
      valores <- as.numeric(fila[!is.na(fila)])  # Excluir valores NA
      if (length(valores) == 0) {
        return(c(NA, NA))  # Devolver NA si no hay datos
      }
      media <- mean(valores, na.rm = TRUE)
      mediana <- median(valores, na.rm = TRUE)
      # Si número par de valores, calcular mediana como media de los dos valores intermedios
      if (length(valores) %% 2 == 0) {
        valores_ordenados <- sort(valores)
        mediana <- mean(valores_ordenados[c(length(valores) / 2, length(valores) / 2 + 1)])
      }
      return(c(media, mediana))
    })
    
    # Convertir a data frame
    resultado_mes <- t(resultado_mes)
    colnames(resultado_mes) <- c("Mean", "Median")
    cbind(Year = datos_mes$Year, resultado_mes)
  })
  
  # Combinar resultados en un solo data frame
  resultado_final <- do.call(cbind, lapply(1:12, function(mes) {
    resultado_mes <- resultados[[mes]]
    colnames(resultado_mes)[-1] <- paste0(colnames(resultado_mes)[-1], "_Month_", mes)
    resultado_mes
  }))
  
  # Eliminar duplicados de la columna Year
  resultado_final <- resultado_final[, !duplicated(colnames(resultado_final))]
  
  # Crear nombres de archivos de salida basados en los datasets usados
  nombre_datasets <- paste0(tools::file_path_sans_ext(basename(archivos_usar)), collapse = "_")
  archivo_median <- file.path(carpeta_inputs, paste0("median_", nombre_datasets, ".txt"))
  archivo_mean <- file.path(carpeta_inputs, paste0("mean_", nombre_datasets, ".txt"))
  
  # Guardar resultados en archivos .txt sin encabezado
  write.table(resultado_final[, c("Year", grep("Mean", colnames(resultado_final), value = TRUE))],
              file = archivo_mean, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  write.table(resultado_final[, c("Year", grep("Median", colnames(resultado_final), value = TRUE))],
              file = archivo_median, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  
  cat("Archivos generados con éxito (sin encabezado):\n")
  cat(" - ", archivo_median, "\n")
  cat(" - ", archivo_mean, "\n")
}

# Llamada a la función
if (dir.exists(CARPETA_INPUTS)) {
  generar_ensemble_mensual(CARPETA_INPUTS, ARCHIVOS_USAR)
} else {
  warning("La carpeta de inputs no existe: ", CARPETA_INPUTS)
}