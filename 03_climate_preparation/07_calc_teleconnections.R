# ==============================================================================
# 12_calc_teleconnections.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Teleconnection Indices Calculator (NAO, WeMO, SOI, MO).
#   Computes daily and monthly atmospheric teleconnection indices from raw surface 
#   pressure (SLP) and temperature series at key nodes.
#
#   Indices Calculated:
#   - NAO (North Atlantic Oscillation): Diff. between Azores/Gibraltar/Lisbon and Reykjavik.
#   - WeMO (Western Mediterranean Oscillation): San Fernando (Cadiz) - Padua (Italy).
#   - SOI (Southern Oscillation Index): Tahiti - Darwin.
#   - MO (Mediterranean Oscillation): Algiers - Cairo.
#
#   Methodology:
#   - Standardizes raw pressure series (Z-scores).
#   - Computes differences between standardised poles (Nodes).
#
#   Inputs:
#   - Directory of Daily Mean Pressure files for each node.
#
#   Outputs:
#   - Daily and Monthly Index Files for each teleconnection pattern.
# ==============================================================================

library(tidyverse)
library(lubridate)

# --- PARÁMETROS A MODIFICAR ---
INPUT_DIR_PRESIONES <- "PLACEHOLDER/path/to/pressure_data"
OUTPUT_DIR_INDICES  <- file.path(INPUT_DIR_PRESIONES, 'telecommunications')
SITIOS_ANALISIS     <- c("argel", "azores", "cairo", "colduzad", "darwin", 
                         "gibraltar", "jaafar", "lisboa", "padoa", "reykjavik", 
                         "sanfer", "tahiti")
# ------------------------------

if (!dir.exists(OUTPUT_DIR_INDICES)) {
  dir.create(OUTPUT_DIR_INDICES, recursive = TRUE)
  message("Carpeta creada: ", OUTPUT_DIR_INDICES)
}

# ================= CARGA Y PREPARACIÓN DE DATOS =================

leer_presion <- function(sitio, ruta_base) {
  archivo <- file.path(ruta_base, paste0(sitio, "_daily_mean.txt"))
  
  if (!file.exists(archivo)) {
    warning(paste("No se encuentra el archivo:", archivo))
    return(NULL)
  }
  
  datos <- read.table(archivo, header = TRUE, stringsAsFactors = FALSE)
  
  datos <- datos %>%
    mutate(date_obj = ymd(as.character(date))) %>%
    mutate(presion_hpa = sfcsi / 100) %>%
    select(date_obj, !!sitio := presion_hpa)
  
  return(datos)
}

message("Leyendo archivos...")
lista_datos <- lapply(SITIOS_ANALISIS, leer_presion, ruta_base = INPUT_DIR_PRESIONES)
lista_datos <- lista_datos[!sapply(lista_datos, is.null)]
df_total <- lista_datos %>% reduce(full_join, by = "date_obj") %>% arrange(date_obj)


message("Datos cargados. Rango: ", min(df_total$date_obj), " a ", max(df_total$date_obj))


# ================= FUNCIONES DE CÁLCULO =================

# Estandarización mensual (Z = (x - mean) / sd)
calcular_anomalias_estandarizadas <- function(df, columnas_sitios) {
  df %>%
    mutate(mes = month(date_obj)) %>%
    group_by(mes) %>%
    mutate(across(all_of(columnas_sitios), 
                  ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE),
                  .names = "{.col}_std")) %>%
    ungroup() %>%
    select(-mes)
}

# Definición de fórmulas (Añade aquí más si necesitas)
calcular_indices <- function(df_std) {
  df_std %>%
    mutate(
      NAO_AzRey  = azores_std - reykjavik_std,
      NAO_GibRey = gibraltar_std - reykjavik_std,
      NAO_LisRey = lisboa_std - reykjavik_std,
      WeMO       = sanfer_std - padoa_std,
      SOI        = tahiti_std - darwin_std,
      MO         = argel_std - cairo_std
    ) %>%
    select(date_obj, NAO_AzRey, NAO_GibRey, NAO_LisRey, WeMO, SOI, MO)
}

# --- FUNCIÓN DE EXPORTACIÓN MODIFICADA ---
exportar_txt <- function(df, nombre_columna, ruta_salida, sufijo) {
  
  # Seleccionamos fecha y el valor específico
  df_export <- df %>% select(date_obj, valor_temp = !!sym(nombre_columna))
  
  # Formateamos fecha a string YYYYMMDD
  df_export$date <- format(df_export$date_obj, "%Y%m%d")
  
  # Filtramos NAs
  df_export <- df_export %>% filter(!is.na(valor_temp))
  
  # Redondeamos a 4 decimales
  df_export$valor_temp <- round(df_export$valor_temp, 4)
  
  # Preparamos la salida renombrando la columna de valor al nombre del índice
  # Esto asegura que la cabecera sea: date [TAB] NombreIndice
  salida <- df_export %>% 
    select(date, !!nombre_columna := valor_temp)
  
  nombre_archivo <- file.path(ruta_salida, paste0(nombre_columna, "_", sufijo, ".txt"))
  
  write.table(salida, file = nombre_archivo, 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# ================= PROCESAMIENTO DIARIO =================
message("Calculando índices diarios...")

df_daily_std <- calcular_anomalias_estandarizadas(df_total, SITIOS_ANALISIS)
indices_diarios <- calcular_indices(df_daily_std)

nombres_indices <- c("NAO_AzRey", "NAO_GibRey", "NAO_LisRey", "WeMO", "SOI", "MO")

for (idx in nombres_indices) {
  exportar_txt(indices_diarios, idx, OUTPUT_DIR_INDICES, "daily")
}

# ================= PROCESAMIENTO MENSUAL =================
message("Calculando índices mensuales...")

# Agregamos presiones (media mensual)
df_monthly <- df_total %>%
  mutate(mes_fecha = floor_date(date_obj, "month")) %>%
  group_by(mes_fecha) %>%
  summarise(across(all_of(SITIOS_ANALISIS), ~ mean(., na.rm = TRUE)), .groups = 'drop') %>%
  rename(date_obj = mes_fecha)

# Estandarizamos sobre la serie mensual
df_monthly_std <- calcular_anomalias_estandarizadas(df_monthly, SITIOS_ANALISIS)
indices_mensuales <- calcular_indices(df_monthly_std)

for (idx in nombres_indices) {
  exportar_txt(indices_mensuales, idx, OUTPUT_DIR_INDICES, "monthly")
}

message("¡Proceso completado! Archivos guardados en: ", OUTPUT_DIR_INDICES)
