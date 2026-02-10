# ==============================================================================
# 10_calc_spei_prec_pet.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Simple SPEI Calculator (From P & PET Files).
#   A streamlined utility to calculate SPEI when both Precipitation (P) and 
#   Potential Evapotranspiration (PET) time series are already available as 
#   separate text files.
#
#   Functionality:
#   - Merges P and PET files by date.
#   - Computes Climatic Water Balance (D = P - PET).
#   - Calculates SPEI for scales 1 to 24 months.
#   - Outputs a comprehensive matrix suitable for drought analysis.
#
#   Inputs:
#   - Precipitation text file.
#   - PET text file.
#
#   Outputs:
#   - SPEI Matrix file (Date + Balance + SPEI_1...SPEI_24).
# ==============================================================================

# --- 1. LIBRERÍAS ---
library(tidyverse)
library(SPEI)

# --- 2. CONFIGURACIÓN Y CARGA DE DATOS ---

path_prec <- 'PLACEHOLDER/path/to/precip_data.txt'
# Asumo que guardaste el PET extraído en esta ruta o similar
path_pet  <- 'PLACEHOLDER/path/to/pet_data.txt'

# Leer archivos
df_prec <- read.table(path_prec, header = TRUE, sep = "\t")
df_pet  <- read.table(path_pet, header = TRUE, sep = "\t")

# --- 3. PREPARACIÓN Y VALIDACIÓN ---

# Unir por fecha para asegurar que cada mes corresponde con su par
# Esto evita errores si falta algún mes en alguno de los dos archivos
df_clima <- inner_join(df_prec, df_pet, by = "date") %>%
  rename(P = precipitation, PET = pet) %>%
  arrange(date)

cat("Datos cargados. N. meses:", nrow(df_clima), "\n")
cat("Periodo:", min(df_clima$date), "-", max(df_clima$date), "\n")

# --- 4. CÁLCULO DEL BALANCE HÍDRICO (D) ---

# D = Precipitation - Potential Evapotranspiration
# Si P < PET, el balance es negativo (déficit de agua)
df_clima$Balance <- df_clima$P - df_clima$PET

# --- 5. CÁLCULO DEL SPEI MULTI-ESCALA (Bucle 1 a 24) ---

# Definir escalas a calcular
max_scale <- 24

# Crear matriz vacía para almacenar resultados
# Filas = Meses, Columnas = Escalas
spei_matrix <- matrix(NA, nrow = nrow(df_clima), ncol = max_scale)
colnames(spei_matrix) <- paste0("SPEI_", 1:max_scale)

cat("Calculando SPEI para escalas 1 a", max_scale, "...\n")

for (i in 1:max_scale) {
  # Calcular SPEI para la escala 'i'
  # Usamos na.rm = TRUE por seguridad, aunque spei() maneja NAs internos
  spei_calc <- spei(df_clima$Balance, scale = i, na.rm = TRUE)
  
  # Extraer la serie numérica (fitted values)
  spei_matrix[, i] <- as.numeric(spei_calc$fitted)
}

# --- 6. EXPORTACIÓN DE RESULTADOS ---

# Unir fechas, balance y la matriz de SPEIs
df_final <- cbind(
  select(df_clima, date, Balance), # Mantenemos fecha y balance original
  as.data.frame(spei_matrix)       # Añadimos las 24 columnas de SPEI
)

# Definir nombre de salida
archivo_salida <- "cdz_spei_completo_1_24.txt"

# Guardar en disco
write.table(df_final, 
            file = archivo_salida, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

cat("==============================================================================\n")
cat("¡Cálculo finalizado!\n")
cat("Archivo guardado como:", archivo_salida, "\n")
cat("Contiene", ncol(df_final), "columnas (Date + Balance + 24 escalas SPEI).\n")
cat("==============================================================================\n")

# --- 7. (OPCIONAL) VISUALIZACIÓN RÁPIDA ---
# Plotear el SPEI de 12 meses para verificar visualmente
plot(df_final$SPEI_12, type="l", col="blue", lwd=1.5,
     main="SPEI escala 12 meses (Verificación)", 
     ylab="SPEI", xlab="Meses desde inicio")
abline(h=0, col="black", lty=2)
abline(h=-1.5, col="red", lty=3) # Umbral de sequía severa