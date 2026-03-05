# ==============================================================================
# 13_plot_monthly_residual.R
# Autor: Marcos Marín-Martín
# Fecha: 2026-02-09
# Descripción: Gráfica de correlaciones entre SPI-12 (Agosto) y RWI.
#              compara múltiples fuentes de datos climáticos.
# ==============================================================================

# Instalar y cargar las librerías necesarias
library(ggplot2)
library(tidyr)
library(showtext)


# --- PARÁMETROS A MODIFICAR ---
OUTPUT_DIR       <- 'PLACEHOLDER/path/to/output_dir'
PLOT_WIDTH       <- 16  # Ancho del gráfico
PLOT_HEIGHT      <- 8  # Alto del gráfico
LINE_WIDTH       <- 1.5  # Anchura de las líneas del gráfico
FONT_FAMILY      <- "Rubik"  # Tipo de letra
LEGEND_TEXT_SIZE <- 18  # Tamaño de la letra de la leyenda

# Rutas a archivos de entrada
INPUT_FILES_LIST <- list(
  TERRA   = "PLACEHOLDER/path/to/res_correlation_terra.txt",
  ROCIO   = "PLACEHOLDER/path/to/res_correlation_rocio.txt",
  SPREAD  = "PLACEHOLDER/path/to/res_correlation_spread.txt",
  EOBS01  = "PLACEHOLDER/path/to/res_correlation_eobs_01.txt",
  EOBS025 = "PLACEHOLDER/path/to/res_correlation_eobs_025.txt",
  ERA5    = "PLACEHOLDER/path/to/res_correlation_era5.txt",
  CRU05   = "PLACEHOLDER/path/to/res_correlation_cru.txt",
  CRU10   = "PLACEHOLDER/path/to/res_correlation_cru_10.txt"
)
# ------------------------------


# Añadir la fuente Rubik
font_add_google(FONT_FAMILY, "rubik")
showtext_auto()

# Crear la carpeta de salida si no existe
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Diccionario completo de colores y etiquetas para los posibles datasets
ALL_COLORS <- c(
  TERRA   = "navajowhite4",
  ROCIO   = "grey13",
  SPREAD  = "yellow1",
  EOBS01  = "turquoise3",
  EOBS025 = "royalblue2",
  ERA5    = "lightpink1",
  CRU05   = "olivedrab2",
  CRU10   = "springgreen3"
)
ALL_LABELS <- c(
  TERRA   = "TerraClimate (obs., res.: 0.04\u00ba)",
  ROCIO   = "ROCIO (AEMET) (obs., res.: 0.045\u00ba)",
  SPREAD  = "SPREAD (obs., res.: 0.045\u00ba)",
  EOBS01  = "E-OBS (obs., res.: 0.1\u00ba)",
  EOBS025 = "E-OBS (obs., res.: 0.25\u00ba)",
  ERA5    = "ERA5 (ren\u00e1lisis, res.: 0.282\u00ba)",
  CRU05   = "CRU (obs., res.: 0.5\u00ba)",
  CRU10   = "CRU (obs., res.: 1.0\u00ba)"
)

# Inicializar una lista vacía para almacenar los datos
data_list <- list()

# Leer cada archivo y extraer la columna de agosto
for (name in names(INPUT_FILES_LIST)) {
  path <- INPUT_FILES_LIST[[name]]

  # Saltar entradas con rutas PLACEHOLDER o que no existen en disco
  if (grepl("PLACEHOLDER", path, fixed = TRUE) || !file.exists(path)) {
    message("Omitiendo dataset '", name, "': archivo no disponible.")
    next
  }

  data <- tryCatch(
    read.table(path, header = TRUE, na.strings = "NA"),
    error = function(e) {
      warning(paste("No se pudo leer el archivo:", path))
      return(NULL)
    }
  )

  if (!is.null(data) && "August" %in% colnames(data)) {
    # Extraer la columna de agosto y limitar a Lag 1 a 24
    august_data <- data$August[1:24]
    # Asegurarse de que tenga exactamente 24 valores, rellenando con NA si es necesario
    if (length(august_data) < 24) {
      august_data <- c(august_data, rep(NA, 24 - length(august_data)))
    }
    data_list[[name]] <- august_data
  } else {
    warning(paste("El archivo", path, "no contiene la columna 'August'. Dataset omitido."))
  }
}

# Verificar que haya al menos un dataset válido
if (length(data_list) == 0) {
  stop("No se encontró ningún archivo de entrada válido. Revisa las rutas en INPUT_FILES_LIST.")
}

valid_names <- names(data_list)

# Crear un data frame dinámico solo con los datasets disponibles
lags <- 1:24
plot_data <- data.frame(Lag = lags, data_list[valid_names])

# Convertir los datos a formato largo para ggplot2
plot_data_long <- pivot_longer(plot_data, cols = -Lag, names_to = "Dataset", values_to = "Value")

# Establecer los niveles del factor en el orden original, solo con los disponibles
plot_data_long$Dataset <- factor(plot_data_long$Dataset, levels = valid_names)

# Filtrar colores y etiquetas a los datasets disponibles
active_colors <- ALL_COLORS[valid_names]
active_labels <- ALL_LABELS[valid_names]

# Crear el gráfico con líneas suaves, dejando los valores NA como huecos
plot <- ggplot(plot_data_long, aes(x = Lag, y = Value, color = Dataset)) +
  geom_line(na.rm = FALSE, linewidth = LINE_WIDTH) +
  scale_color_manual(values = active_colors, labels = active_labels) +
  labs(title = "Pearson correlation between SPI-12 and RWI for August",
       x = "Lag (months)",
       y = "R",
       color = "Dataset") +
  theme_minimal() +
  theme(text = element_text(family = "rubik", color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),
        legend.text = element_text(size = LEGEND_TEXT_SIZE, color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title = element_text(color = "black"))

# Guardar el gráfico como PDF y PNG con dimensiones 2:1
ggsave(filename = file.path(OUTPUT_DIR, "correlation_plot.pdf"), plot = plot, width = PLOT_WIDTH, height = PLOT_HEIGHT)
ggsave(filename = file.path(OUTPUT_DIR, "correlation_plot.png"), plot = plot, width = PLOT_WIDTH, height = PLOT_HEIGHT)

