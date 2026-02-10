# ==============================================================================
# 13_sweep_seasonal_nao.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Seasonal NAO Signal Sweeper.
#   Performs an exhaustive correlation analysis ("Sweep") between a target 
#   reconstruction variable (e.g., precipitation, tree-ring width) and the NAO 
#   index across multiple seasonal windows.
#
#   Goal:
#   - To identify the specific seasonal window (e.g., Winter DJF, Extended Winter 
#     ONDJFM, Spring MAM) where the NAO signal is strongest in the proxy record.
#
#   Visualisation:
#   - Generates a Heatmap of correlations (Variable vs. Reconstruction Site).
#   - Highlights significant correlations (* p<0.05, ** p<0.01).
#
#   Inputs:
#   - Reconstruction/Precipitation file.
#   - Monthly NAO index file.
#
#   Outputs:
#   - Correlation Matrix Heatmap.
#   - Console summary of best performing windows.
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

# --- PARÁMETROS A MODIFICAR ---
RUTA_PRECIP <- "PLACEHOLDER/path/to/precipitation.txt"
RUTA_NAO    <- "PLACEHOLDER/path/to/nao_index.txt"
# ------------------------------

# --- 1. CARGA DE DATOS ---
# Cargar NAO
nao_raw <- read.table(RUTA_NAO, header = FALSE, fill = TRUE)
colnames(nao_raw) <- c("Year", month.abb) # Year, Jan, Feb...

# Cargar Precipitación
precip <- read.table(RUTA_PRECIP, header = TRUE, sep = "\t") # Ajusta sep="" si es necesario

# --- 2. INGENIERÍA DE VARIABLES (CREAR ESTACIONES) ---

# Necesitamos el Diciembre del año anterior para calcular el invierno (DJF)
# Creamos una columna 'Dec_prev' desplazando la columna Dec una posición abajo
nao_raw$Dec_prev <- lag(nao_raw$Dec, 1)

# Calculamos medias estacionales (puedes añadir más combinaciones aquí)
nao_calc <- nao_raw %>%
  rowwise() %>%
  mutate(
    # Estaciones Clásicas
    Win_DJF = mean(c(Dec_prev, Jan, Feb), na.rm = TRUE), # Invierno Meteorológico
    Win_JFM = mean(c(Jan, Feb, Mar), na.rm = TRUE),      # Invierno Tardío
    Spr_MAM = mean(c(Mar, Apr, May), na.rm = TRUE),      # Primavera
    
    # Ventanas largas
    Wet_ONDJFM = mean(c(Oct, Nov, Dec, Jan, Feb, Mar), na.rm = TRUE), # Semestre húmedo
    
    # Meses individuales clave (usamos los ya existentes)
  ) %>%
  select(Year, Jan, Feb, Mar, Apr, May, Jun, Oct, Nov, Dec, 
         Win_DJF, Win_JFM, Spr_MAM, Wet_ONDJFM) %>%
  ungroup()

# Unir con precipitación
df_total <- inner_join(precip, nao_calc, by = "Year") %>% na.omit()

# --- 3. CÁLCULO DE CORRELACIONES (BARRIDO) ---

# Definimos las variables de NAO a testear
vars_nao <- c("Jan", "Feb", "Mar", "Apr", "May", "Oct", "Nov", "Dec", 
              "Win_DJF", "Win_JFM", "Spr_MAM", "Wet_ONDJFM")

# Función para calcular correlación y p-valor
calcular_correlaciones <- function(serie_precip, nombre_sitio) {
  resultados <- data.frame(Variable = character(), Cor = numeric(), P_val = numeric(), Sitio = character())
  
  for (v in vars_nao) {
    test <- cor.test(df_total[[v]], serie_precip)
    resultados <- rbind(resultados, data.frame(
      Variable = v,
      Cor = test$estimate,
      P_val = test$p.value,
      Sitio = nombre_sitio
    ))
  }
  return(resultados)
}

cor_atlas  <- calcular_correlaciones(df_total$atlas, "Atlas")
cor_teruel <- calcular_correlaciones(df_total$teruel, "Teruel")

df_corr <- rbind(cor_atlas, cor_teruel)

# Crear columna de significancia para el gráfico (* si p < 0.05)
df_corr$Signif <- ifelse(df_corr$P_val < 0.05, "*", "")
df_corr$Signif_Strong <- ifelse(df_corr$P_val < 0.01, "**", df_corr$Signif)

# Ordenar los factores para que salgan en orden lógico en el gráfico
df_corr$Variable <- factor(df_corr$Variable, levels = vars_nao)

# --- 4. VISUALIZACIÓN: HEATMAP DE CORRELACIONES ---

p_heat <- ggplot(df_corr, aes(x = Variable, y = Sitio, fill = Cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f%s", Cor, Signif)), color = "black", size = 4) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, limit = c(-0.6, 0.6)) +
  labs(title = "Matriz de Sensibilidad a la NAO (1822-2024)",
       subtitle = "Rojo = Correlación Negativa (NAO+ -> Seco). * = Signif p<0.05",
       x = "Ventana Temporal de la NAO", y = "Reconstrucción", fill = "Corr (r)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_heat)

# --- 5. RESULTADOS EN TEXTO ---
cat("\n--- MEJORES CORRELACIONES ENCONTRADAS ---\n")
top_atlas <- df_corr %>% filter(Sitio == "Atlas") %>% arrange(P_val) %>% head(1)
top_teruel <- df_corr %>% filter(Sitio == "Teruel") %>% arrange(P_val) %>% head(1)

cat(sprintf("Mejor señal para ATLAS: %s (r = %.3f, p = %.4f)\n", top_atlas$Variable, top_atlas$Cor, top_atlas$P_val))
cat(sprintf("Mejor señal para TERUEL: %s (r = %.3f, p = %.4f)\n", top_teruel$Variable, top_teruel$Cor, top_teruel$P_val))
