# ==============================================================================
# SCRIPT MAESTRO: DENDROTOOLS HÍBRIDO (SELECTOR RES/STD + STATS + PLOTS VALIDACIÓN)
# AUTOR:    Marcos Marín-Martín
# FECHA:    Febrero 2026 (Actualizado para nueva estructura)
#
# DESCRIPCIÓN:
#   1. Detecta automáticamente si el clima es DIARIO o MENSUAL.
#   2. Permite ELEGIR entre cronología RESIDUAL o ESTÁNDAR.
#   3. Ejecuta daily_response() o monthly_response().
#   4. Devuelve y GUARDA los estadísticos de validación (RE, CE).
#   5. GENERA y GUARDA el gráfico de Calibración-Verificación (Split-sample).
# ==============================================================================

# 1. PARÁMETROS A MODIFICAR
# -------------------------------------------------------------------
# Rutas ajustadas a la nueva ubicación: PLACEHOLDER/path/to/dendroclima
ruta_datos_clima <- 'PLACEHOLDER/path/to/climate_data.txt'
ruta_cronologia  <- 'PLACEHOLDER/path/to/chronology.txt'
ruta_salida      <- 'PLACEHOLDER/path/to/output_correlations/'
alpha            <- 0.05

# --- SELECCIÓN DE TIPO DE CRONOLOGÍA ---
# Opciones: "res" (Residual) o "std" (Estándar)
TIPO_CRONO <- "res"  

# --- FILTRO DE AÑOS ---
FILTRAR_ANIOS <- FALSE
ANIO_INICIO   <- 1936
ANIO_FIN      <- 2019

# 2. GENERACIÓN DE IDENTIFICADORES
# -------------------------------------------------------------------
nombre_crono_base <- tools::file_path_sans_ext(basename(ruta_cronologia))
nombre_clima_base <- tools::file_path_sans_ext(basename(ruta_datos_clima))

# Añadimos el tipo de crono al nombre para no mezclar archivos
id_analisis <- paste0(nombre_crono_base, "_", TIPO_CRONO, "_vs_", nombre_clima_base)

if (FILTRAR_ANIOS) {
  id_analisis <- paste0(id_analisis, "_", ANIO_INICIO, "-", ANIO_FIN)
}

message(">>> ID del análisis: ", id_analisis)
message(">>> Tipo de serie seleccionada: ", toupper(TIPO_CRONO))

# 3. CARGA DE LIBRERÍAS
# -------------------------------------------------------------------
library(dendroTools)
library(dplR)
library(ggplot2)
library(dplyr)
library(reshape2)
library(colourvalues)
library(gridExtra)
library(lmtest)

# Forzar que select y otros vengan de dplyr en los pipes
select <- dplyr::select

dirs <- c("txt", "plots", "matrices")
for(d in dirs) dir.create(file.path(ruta_salida, d), recursive = TRUE, showWarnings = FALSE)

# 4. PROCESAMIENTO INTELIGENTE DEL CLIMA (DETECCIÓN DE PERIODICIDAD)
# -------------------------------------------------------------------
datos_clima <- read.table(ruta_datos_clima, header = TRUE, sep = "\t")
colnames(datos_clima) <- gsub('"', '', colnames(datos_clima))

# Detectar columna de valor (la que no es 'date')
vars_cols <- names(datos_clima)
col_variable <- vars_cols[!vars_cols %in% c("date")][1]
message(">>> Variable climática detectada: ", col_variable)
datos_clima$valor <- datos_clima[[col_variable]]

# --- ALGORITMO DE DETECCIÓN MEJORADO ---
sample_date_str <- as.character(datos_clima$date[1])
n_char <- nchar(sample_date_str)

TIPO_DATOS <- "DESCONOCIDO"

if (n_char == 6) {
  TIPO_DATOS <- "MENSUAL"
  message(">>> DETECTADO: Datos MENSUALES (Por formato YYYYMM)")
} else if (n_char == 8) {
  check_n <- min(nrow(datos_clima), 20)
  fechas_check <- as.Date(as.character(datos_clima$date[1:check_n]), format="%Y%m%d")
  diff_days <- median(as.numeric(diff(fechas_check)), na.rm = TRUE)
  
  if (diff_days > 25) {
    TIPO_DATOS <- "MENSUAL"
    message(">>> DETECTADO: Datos MENSUALES (Por periodicidad)")
  } else {
    TIPO_DATOS <- "DIARIO"
    message(">>> DETECTADO: Datos DIARIOS (Por periodicidad)")
  }
} else {
  stop("ERROR: Formato de fecha desconocido (ni 6 ni 8 caracteres).")
}

# --- CONSTRUCCIÓN DE LA MATRIZ ---
if (TIPO_DATOS == "MENSUAL") {
  datos_clima$year  <- as.numeric(substr(as.character(datos_clima$date), 1, 4))
  datos_clima$month <- as.numeric(substr(as.character(datos_clima$date), 5, 6))
  env_data <- dcast(datos_clima, year ~ month, value.var = "valor")
  for (m in 1:12) {
    if (!as.character(m) %in% colnames(env_data)) env_data[[as.character(m)]] <- NA
  }
  env_data <- env_data[, c("year", as.character(1:12))]
  
} else {
  datos_clima$date_fmt <- as.Date(as.character(datos_clima$date), format = "%Y%m%d")
  datos_clima$year <- format(datos_clima$date_fmt, "%Y")
  datos_clima$day_of_year <- as.numeric(format(datos_clima$date_fmt, "%j"))
  env_data <- dcast(datos_clima, year ~ day_of_year, value.var = "valor")
  max_days <- 366
  if (ncol(env_data) < max_days + 1) {
    env_data <- cbind(env_data, matrix(NA, nrow = nrow(env_data), ncol = max_days - ncol(env_data) + 1))
  }
  colnames(env_data) <- c("year", paste0("X", 1:max_days))
}

# 5. PROCESAMIENTO CRONOLOGÍA (CON SELECTOR)
# -------------------------------------------------------------------
response <- read.table(ruta_cronologia, header = TRUE, sep = "\t")

# Normalizar nombre del año
if(!"Year" %in% names(response)) names(response)[1] <- "Year"

# SELECCIÓN ESTRICTA SEGÚN INTERRUPTOR
if (TIPO_CRONO == "res") {
  if (!"res" %in% names(response)) stop("ERROR CRÍTICO: Has seleccionado 'res' pero el archivo no tiene esa columna.")
  response <- response %>% dplyr::select(Year, res)
  names(response)[2] <- "TRW"
  
} else if (TIPO_CRONO == "std") {
  if (!"std" %in% names(response)) stop("ERROR CRÍTICO: Has seleccionado 'std' pero el archivo no tiene esa columna.")
  response <- response %>% dplyr::select(Year, std)
  names(response)[2] <- "TRW"
  
} else {
  stop("ERROR: El parámetro TIPO_CRONO debe ser 'res' o 'std'.")
}

response$Year <- as.numeric(response$Year)

# ==============================================================================
# 6. ALINEACIÓN, FILTRADO Y ACTUALIZACIÓN DEL ID
# ==============================================================================
common_years <- intersect(response$Year, env_data$year)

# Aplicar filtro manual si está activado
if (FILTRAR_ANIOS) {
  target_years <- ANIO_INICIO:ANIO_FIN
  common_years <- intersect(common_years, target_years)
  if (length(common_years) == 0) stop("ERROR: Rango de años seleccionado sin datos comunes.")
}

if(length(common_years) < 10) stop("Error: Demasiados pocos años en común (<10).")

# --- ACTUALIZACIÓN DINÁMICA DEL NOMBRE DEL ARCHIVO ---
rango_real_inicio <- min(common_years)
rango_real_fin    <- max(common_years)
etiqueta_anios    <- paste0(rango_real_inicio, "-", rango_real_fin)

id_analisis <- paste0(nombre_crono_base, "_", TIPO_CRONO, "_vs_", nombre_clima_base, "_", etiqueta_anios)

message(">>> Periodo efectivo de análisis: ", etiqueta_anios)
message(">>> ID Actualizado para archivos: ", id_analisis)

# Recortar datos
response_final <- response[response$Year %in% common_years, ]
env_data_final <- env_data[env_data$year %in% common_years, ]

row.names(response_final) <- response_final$Year
response_final$Year <- NULL

row.names(env_data_final) <- env_data_final$year
env_data_final$year <- NULL

# 7. EJECUCIÓN DEL ANÁLISIS
# -------------------------------------------------------------------
message(">>> Calculando correlaciones (Modo ", TIPO_DATOS, ")...")

if (TIPO_DATOS == "DIARIO") {
  resultados <- daily_response(
    alpha = alpha, 
    response = response_final, 
    env_data = env_data_final,
    method = "cor", metric = "r.squared", cor_method = "pearson", 
    lower_limit = 1, upper_limit = 732, 
    row_names_subset = TRUE, remove_insignificant = TRUE, 
    previous_year = TRUE, reference_window = "end", aggregate_function = "sum"
  )
} else {
  resultados <- monthly_response(
    alpha = alpha,
    response = response_final,
    env_data = env_data_final,
    method = "cor", metric = "r.squared", cor_method = "pearson",
    lower_limit = 1, upper_limit = 24, 
    row_names_subset = TRUE, remove_insignificant = TRUE,
    previous_year = TRUE, reference_window = "end", aggregate_function = "sum"
  )
}

# ==============================================================================
# 8. CÁLCULOS DE VENTANA ÓPTIMA Y ESTADÍSTICOS SPLIT-SAMPLE
# ==============================================================================

colors <- colourvalues::colour_values(1:100, palette = "diverge_hsv")
resultados$plot_heatmap <- resultados$plot_heatmap +
  scale_fill_gradientn(colors = colors) +
  ggtitle(paste0(id_analisis, " (", TIPO_DATOS, ")"))

alpha_str <- sprintf("alfa%02d", alpha * 100)
sufijo_archivo <- paste0(id_analisis, "_", alpha_str)

message(">>> Guardando Heatmap...")
ggsave(
  filename = file.path(ruta_salida, "plots", paste0("Heatmap_", sufijo_archivo, ".pdf")),
  plot = resultados$plot_heatmap,
  width = 10, 
  height = 8
)

# --- A. DETECTAR VENTANA ÓPTIMA (PARA EL ENCABEZADO) ---
matriz_stats <- resultados$calculations
max_valor <- max(matriz_stats, na.rm = TRUE)
pos_max <- which(matriz_stats == max_valor, arr.ind = TRUE)[1, ]
duracion_optima <- as.numeric(row.names(matriz_stats)[pos_max[1]])
indice_inicio   <- pos_max[2]
nombre_inicio   <- colnames(matriz_stats)[indice_inicio]

unidad_tiempo <- ifelse(TIPO_DATOS == "MENSUAL", "meses", "días")
texto_ventana_optima <- paste0("Inicio: ", nombre_inicio, " | Duración: ", duracion_optima, " ", unidad_tiempo)
texto_estadistico    <- paste0("Mejor R global: ", round(max_valor, 4))

write.table(as.data.frame(resultados$calculations), 
            file.path(ruta_salida, "txt", paste0("Correl_Matrix_", sufijo_archivo, ".txt")), 
            sep = "\t", row.names = TRUE, col.names = NA)

# ==============================================================================
# 9. VALIDACIÓN CRUZADA (SPLIT-SAMPLE) Y GENERACIÓN DE REPORTE ÚNICO
# ==============================================================================

tabla_split <- NULL
tabla_dendro <- resultados$cross_validation

if (!is.null(resultados$optimised_return)) {
  
  message(">>> Calculando estadísticas Split-Sample...")
  
  datos_modelo <- resultados$optimised_return
  names(datos_modelo) <- c("TRW", "Clima")
  datos_modelo$Year <- as.numeric(row.names(datos_modelo))
  
  mitad <- floor((min(datos_modelo$Year) + max(datos_modelo$Year)) / 2)
  early_period <- subset(datos_modelo, Year <= mitad)
  late_period  <- subset(datos_modelo, Year > mitad)
  
  calc_dendro_stats <- function(obs, pred, cal_mean) {
    rmse <- sqrt(mean((obs - pred)^2))
    RE <- 1 - (sum((obs - pred)^2) / sum((obs - cal_mean)^2))
    CE <- 1 - (sum((obs - pred)^2) / sum((obs - mean(obs))^2))
    return(c(RMSE = rmse, RE = RE, CE = CE))
  }
  
  mod_early <- lm(Clima ~ TRW, data = early_period)
  pred_late <- predict(mod_early, newdata = late_period)
  skills_late <- calc_dendro_stats(late_period$Clima, pred_late, mean(early_period$Clima))
  
  mod_late <- lm(Clima ~ TRW, data = late_period)
  pred_early <- predict(mod_late, newdata = early_period)
  skills_early <- calc_dendro_stats(early_period$Clima, pred_early, mean(late_period$Clima))
  
  dw_early <- dwtest(mod_early)
  dw_late  <- dwtest(mod_late)
  
  tabla_split <- data.frame(
    Model = c("Calib_Early_Ver_Late", "Calib_Late_Ver_Early"),
    Period_Cal = c(paste(min(early_period$Year), max(early_period$Year), sep="-"),
                   paste(min(late_period$Year), max(late_period$Year), sep="-")),
    Period_Ver = c(paste(min(late_period$Year), max(late_period$Year), sep="-"),
                   paste(min(early_period$Year), max(early_period$Year), sep="-")),
    R2_Cal = round(c(summary(mod_early)$r.squared, summary(mod_late)$r.squared), 3),
    R2_Ver = round(c(cor(late_period$Clima, pred_late)^2, cor(early_period$Clima, pred_early)^2), 3),
    RE = round(c(skills_late["RE"], skills_early["RE"]), 3),
    CE = round(c(skills_late["CE"], skills_early["CE"]), 3),
    DW_Stat = round(c(dw_early$statistic, dw_late$statistic), 3),
    DW_pVal = round(c(dw_early$p.value, dw_late$p.value), 3)
  )
  
  df_A <- data.frame(Year = early_period$Year, Obs = early_period$Clima, Est = pred_early)
  df_B <- data.frame(Year = late_period$Year, Obs = late_period$Clima, Est = pred_late)
  
  plot_split_fn <- function(df, tit, r2) {
    ggplot(df, aes(x = Year)) +
      geom_line(aes(y = Obs, color = "Obs"), linewidth = 0.6) +
      geom_line(aes(y = Est, color = "Est"), linewidth = 0.6, linetype = "dashed") +
      annotate("text", x = mean(df$Year), y = max(df$Obs), label = paste0("r2=", round(r2,2)), fontface="bold") +
      scale_color_manual(name="", values=c("Obs"="black", "Est"="red")) +
      labs(title = tit, y = col_variable) + theme_bw() + theme(legend.position="bottom")
  }
  
  p1 <- plot_split_fn(df_A, "Verificación: Periodo Antiguo", cor(df_A$Obs, df_A$Est)^2)
  p2 <- plot_split_fn(df_B, "Verificación: Periodo Reciente", cor(df_B$Obs, df_B$Est)^2)
  
  pdf(file.path(ruta_salida, "plots", paste0("Split_Verification_", sufijo_archivo, ".pdf")), width=10, height=5)
  grid.arrange(p1, p2, ncol=2)
  dev.off()
}

# ==============================================================================
# 10. ESCRITURA DEL ARCHIVO UNIFICADO FINAL
# ==============================================================================

archivo_reporte <- file.path(ruta_salida, "txt", paste0("Validation_Stats_", sufijo_archivo, ".txt"))
con <- file(archivo_reporte, "w")

writeLines("======================================================", con)
writeLines(paste("ANÁLISIS:", id_analisis), con)
writeLines(paste("FECHA:   ", Sys.time()), con)
writeLines("======================================================", con)
writeLines("RESUMEN DE VENTANA CLIMÁTICA ÓPTIMA:", con)
writeLines(paste(" >>", texto_estadistico), con)
writeLines(paste(" >>", texto_ventana_optima), con)
if (TIPO_DATOS == "MENSUAL") {
  writeLines("    (Nota: Si indice > 12, corresponde a meses del año actual)", con)
}
writeLines("======================================================\n", con)

if (!is.null(tabla_split)) {
  writeLines("--- SECCIÓN 1: VALIDACIÓN CRUZADA (SPLIT-SAMPLE) ---", con)
  writeLines("    (Calculado dividiendo la serie en dos mitades)", con)
  write.table(tabla_split, con, sep = "\t", row.names = FALSE, quote = FALSE)
  writeLines("\n", con)
}

if (!is.null(tabla_dendro)) {
  writeLines("--- SECCIÓN 2: ESTADÍSTICOS DENDROTOOLS (BOOTSTRAP/K-FOLD) ---", con)
  writeLines("    (Generados automáticamente por la librería)", con)
  write.table(tabla_dendro, con, sep = "\t", row.names = FALSE, quote = FALSE)
}

close(con)

message("\n>>> REPORTE GENERADO EN: ", archivo_reporte)
cat("\n--- RESUMEN ---\n")
cat(texto_ventana_optima, "\n")
cat(texto_estadistico, "\n")
if(!is.null(tabla_split)) print(tabla_split)
