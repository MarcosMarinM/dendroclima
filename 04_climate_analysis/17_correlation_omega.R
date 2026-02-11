# ==============================================================================
# 17_correlation_omega.R
# Author: Marcos Marín-Martín
# Date: 2026-02-11
# Description:
#   Optimisation of Omega (Vertical Velocity) Correlations.
#   Iterates through pressure levels and seasonal windows to find the optimal
#   configuration where Omega best correlates (typically negatively) with tree growth.
# ==============================================================================

library(ncdf4); library(dplyr); library(ggplot2); library(reshape2)

# --- 1. CONFIGURACIÓN ---
# RUTA AL ARCHIVO OMEGA
file_omega <- "PLACEHOLDER/path/to/omega.mon.mean.nc"
# RUTA A LOS ÁRBOLES
file_tree  <- "PLACEHOLDER/path/to/chronology.txt"

# CAJA ATLAS (Misma que usamos para el viento)
lat_box <- c(30, 37); lon_box <- c(350, 5)

# TEMPORADAS A TESTEAR
seasons <- list("NDJ"=c(11,12,1), "DJF"=c(12,1,2), "JFM"=c(1,2,3), 
                "FMA"=c(2,3,4), "MAM"=c(3,4,5), "AMJ"=c(4,5,6))

# --- 2. CARGA DE DATOS ---
print("Cargando NetCDF Omega...")
nc <- nc_open(file_omega)
levels <- ncvar_get(nc, "level")
lat <- ncvar_get(nc, "lat")
lon <- ncvar_get(nc, "lon")
time <- ncvar_get(nc, "time")
dates <- data.frame(date = as.Date(time/24, origin="1800-01-01"))
dates$year <- as.numeric(format(dates$date, "%Y"))
dates$month <- as.numeric(format(dates$date, "%m"))

# Indices espaciales
idx_lat <- which(lat >= lat_box[1] & lat <= lat_box[2])
idx_lon <- c(which(lon >= lon_box[1]), which(lon <= lon_box[2]))

# Matriz para guardar resultados
results_mat <- matrix(NA, nrow=length(levels), ncol=length(seasons))
rownames(results_mat) <- levels
colnames(results_mat) <- names(seasons)

tree_data <- read.table(file_tree, header = TRUE)

# --- 3. BUCLE DE OPTIMIZACIÓN (NIVEL x TEMPORADA) ---
print("Analizando todos los niveles y temporadas...")

for(i in 1:length(levels)) {
  # Extraer nivel i
  w_slice <- ncvar_get(nc, "omega", start=c(1,1,i,1), count=c(-1,-1,1,-1))
  # Recortar Atlas
  w_atlas <- w_slice[idx_lon, idx_lat, ]
  # Serie mensual promedio
  dates$val <- apply(w_atlas, 3, mean, na.rm=TRUE)
  
  for(j in 1:length(seasons)) {
    s_m <- seasons[[j]]
    # Ajuste año climático
    df_temp <- dates
    if(any(s_m %in% c(10,11,12))) {
      df_temp$clim_year <- ifelse(df_temp$month %in% c(10,11,12), df_temp$year+1, df_temp$year)
    } else {
      df_temp$clim_year <- df_temp$year
    }
    
    # Promedio estacional
    seas_avg <- df_temp %>% 
      filter(month %in% s_m) %>% 
      group_by(clim_year) %>% 
      summarise(omega_mean=mean(val), n=n()) %>% 
      filter(n==3)
    
    # Correlación
    merged <- inner_join(tree_data, seas_avg, by=c("year"="clim_year"))
    if(nrow(merged)>20) {
      # Guardamos la correlación
      results_mat[i,j] <- cor(merged$res, merged$omega_mean)
    }
  }
}
nc_close(nc)

# --- 4. RESULTADOS: EL MAPA DE CALOR ---
melted <- melt(results_mat)
colnames(melted) <- c("Level", "Season", "r")
melted$Level <- factor(melted$Level, levels=sort(levels, decreasing=TRUE)) # Ordenar visualmente

# Encontrar el MEJOR (El más negativo es el mejor para Omega)
best_idx <- which.min(melted$r) # Buscamos el mínimo (correlación negativa más fuerte)
best_combo <- melted[best_idx, ]

# Gráfico Matriz
p1 <- ggplot(melted, aes(x=Season, y=Level, fill=r)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       name="Pearson r") +
  geom_text(aes(label=round(r, 2)), size=2.5, color="black") +
  # Marcar el ganador con un cuadro amarillo
  geom_tile(data=best_combo, aes(x=Season, y=Level), fill=NA, color="yellow", size=1.5) +
  labs(title="Matriz Omega: Buscando el AZUL intenso (Lluvia)", 
       subtitle=paste("Ganador:", best_combo$Season, "a", best_combo$Level, "hPa | r =", round(best_combo$r, 3))) +
  theme_minimal()

print(p1)

# --- 5. GRÁFICO DE RECONSTRUCCIÓN DEL GANADOR ---
cat("\nGANADOR ABSOLUTO:\n")
cat(paste("Nivel:", best_combo$Level, "hPa\n"))
cat(paste("Temporada:", best_combo$Season, "\n"))
cat(paste("Correlación:", round(best_combo$r, 4), "(Negativa es buena)\n"))

# Recalcular serie del ganador para plotear
nc <- nc_open(file_omega)
lev_idx <- which(abs(levels - as.numeric(as.character(best_combo$Level))) < 1)
w_slice <- ncvar_get(nc, "omega", start=c(1,1,lev_idx,1), count=c(-1,-1,1,-1))
nc_close(nc)
w_atlas <- w_slice[idx_lon, idx_lat, ]
dates$val <- apply(w_atlas, 3, mean, na.rm=TRUE)

s_m <- seasons[[as.character(best_combo$Season)]]
df_temp <- dates
if(any(s_m %in% c(10,11,12))) {
  df_temp$clim_year <- ifelse(df_temp$month %in% c(10,11,12), df_temp$year+1, df_temp$year)
} else { df_temp$clim_year <- df_temp$year }

final_series <- df_temp %>% filter(month %in% s_m) %>% group_by(clim_year) %>% 
  summarise(omega_mean=mean(val))

merged_final <- inner_join(tree_data, final_series, by=c("year"="clim_year"))

# Modelo lineal (Omega predice Arboles)
model <- lm(res ~ omega_mean, data=merged_final)
r_sq <- summary(model)$r.squared

# Graficar (INVERTIMOS OMEGA VISUALMENTE PARA QUE SE VEA LA SINCRONÍA)
ggplot(merged_final, aes(x=year)) +
  # Omega invertido y escalado
  geom_line(aes(y=scale(omega_mean)*-1, color="Omega (Invertido: Arriba=Lluvia)"), size=1) +
  # Árboles escalados
  geom_line(aes(y=scale(res), color="Crecimiento Árboles"), alpha=0.7) +
  scale_color_manual(values=c("black", "blue")) +
  labs(title=paste("Reconstrucción: Árboles vs Omega (", best_combo$Season, best_combo$Level, "hPa)"),
       subtitle=paste("Correlación r =", round(best_combo$r, 3), "| Varianza Explicada R2 =", round(r_sq*100,1), "%"),
       y="Anomalías Estandarizadas (Z-score)") +
  theme_minimal() + theme(legend.position="bottom")
