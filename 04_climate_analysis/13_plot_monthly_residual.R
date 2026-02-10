# ==============================================================================
# 02_grafico_mes_cronoresidual.R
# Autor: Marcos Marín-Martín
# Fecha: 2026-02-09
# Descripción: Gráfica de correlaciones entre SPI-12 (Agosto) y RWI.
#              compara múltiples fuentes de datos climáticos.
# ==============================================================================

# Instalar y cargar las librerías necesarias
# Instalar y cargar las librerías necesarias
library(ggplot2)
library(tidyr)
library(showtext)


# --- PARÁMETROS A MODIFICAR ---
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

# Inicializar una lista vacía para almacenar los datos
data_list <- list()

# Leer cada archivo y extraer la columna de agosto
for (name in names(INPUT_FILES_LIST)) {
  data <- tryCatch(
    read.table(INPUT_FILES_LIST[[name]], header = TRUE, na.strings = "NA"),
    error = function(e) {
      warning(paste("No se pudo leer el archivo:", INPUT_FILES_LIST[[name]]))
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
  } else {
    warning(paste("El archivo", INPUT_FILES_LIST[[name]], "no contiene la columna 'August'."))
    august_data <- rep(NA, 24)
  }
  
  data_list[[name]] <- august_data
}

# Crear un data frame para graficar
lags <- 1:24
plot_data <- data.frame(
  Lag = lags,
  TERRA = data_list$TERRA,
  ROCIO = data_list$ROCIO,
  SPREAD = data_list$SPREAD,
  EOBS01 = data_list$EOBS01,
  EOBS025 = data_list$EOBS025,
  ERA5 = data_list$ERA5,
  CRU05 = data_list$CRU05,
  CRU10 = data_list$CRU10
)

# Convertir los datos a formato largo para ggplot2
plot_data_long <- pivot_longer(plot_data, cols = -Lag, names_to = "Dataset", values_to = "Value")

# Establecer los niveles del factor para el orden deseado en la leyenda
plot_data_long$Dataset <- factor(plot_data_long$Dataset, levels = c(
  "TERRA", "ROCIO", "SPREAD", "EOBS01", "EOBS025", "ERA5", "CRU05", "CRU10"
))

# Crear el gráfico con líneas suaves, dejando los valores NA como huecos
plot <- ggplot(plot_data_long, aes(x = Lag, y = Value, color = Dataset)) +
  geom_line(na.rm = FALSE, linewidth = LINE_WIDTH) +
  scale_color_manual(values = c(
    "TERRA" = "navajowhite4",
    "ROCIO" = "grey13",
    "SPREAD" = "yellow1",
    "EOBS01" = "turquoise3",
    "EOBS025" = "royalblue2",
    "ERA5" = "lightpink1",
    "CRU05" = "olivedrab2",
    "CRU10" = "springgreen3"
  ),
  labels = c(
    "TerraClimate (obs., res.: 0.04º)",
    "ROCIO (AEMET) (obs., res.: 0.045º)",
    "SPREAD (obs., res.: 0.045º)",
    "E-OBS (obs., res.: 0.1º)",
    "E-OBS (obs., res.: 0.25º)",
    "ERA5 (reanálisis, res.: 0.282º)",
    "CRU (obs., res.: 0.5º)",
    "CRU (obs., res.: 1.0º)"
  )) +
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