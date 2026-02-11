# ==============================================================================
# 15_field_correlations.R
# Author: Marcos Marín-Martín
# Date: 2026-02-11
# Description:
#   Field Correlations Analysis (Wind at 300hPa).
#   Calculates field correlations between tree-ring chronologies and zonal wind,
#   identifying dipole centers (Green/Yellow crosses) of maximum correlation.
# ==============================================================================

library(ncdf4)
library(dplyr)
library(ggplot2)
library(maps)

# --- CONFIGURACIÓN ---
file_nc   <- "PLACEHOLDER/path/to/uwnd.mon.mean.nc"
file_tree <- "PLACEHOLDER/path/to/chronology.txt"
nivel_hpa <- 300  # Nivel exploratorio inicial

# --- CARGA Y PROCESAMIENTO ---
print("Cargando y procesando datos...")
nc <- nc_open(file_nc)
lev_idx <- which(abs(ncvar_get(nc, "level") - nivel_hpa) < 1)
uwnd <- ncvar_get(nc, "uwnd", start=c(1,1,lev_idx,1), count=c(-1,-1,1,-1))
lat <- ncvar_get(nc, "lat"); lon <- ncvar_get(nc, "lon")
dates <- as.Date(ncvar_get(nc, "time")/24, origin="1800-01-01")
nc_close(nc)

# Calcular Invierno (DJF)
years <- as.numeric(format(dates, "%Y"))
months <- as.numeric(format(dates, "%m"))
unique_years <- unique(years)[-1] 

winter_wind <- array(NA, dim=c(length(lon), length(lat), length(unique_years)))

for(i in 1:length(unique_years)) {
  yr <- unique_years[i]
  # Dic año previo + Ene + Feb año actual
  idx <- which((years == (yr-1) & months == 12) | (years == yr & months %in% 1:2))
  if(length(idx) == 3) winter_wind[,,i] <- apply(uwnd[,,idx], c(1,2), mean)
}

# --- CORRELACIÓN ---
print("Calculando matriz de correlaciones...")
tree_data <- read.table(file_tree, header = TRUE)
comm_yrs <- intersect(tree_data$year, unique_years)

# Indices coincidentes
idx_t <- match(comm_yrs, tree_data$year)
idx_w <- match(comm_yrs, unique_years)

cor_map <- matrix(NA, nrow=length(lon), ncol=length(lat))
tree_vec <- tree_data$res[idx_t]

for(i in 1:length(lon)) {
  for(j in 1:length(lat)) {
    if(!any(is.na(winter_wind[i,j,idx_w]))) {
      cor_map[i,j] <- cor(tree_vec, winter_wind[i,j,idx_w])
    }
  }
}

# --- DETECTAR POLOS (CRUCES) ---
# Encontrar coordenadas del Máximo Positivo (Cruz Verde)
max_val <- max(cor_map, na.rm=TRUE)
idx_max <- which(cor_map == max_val, arr.ind = TRUE)
lon_max <- lon[idx_max[1]]
lat_max <- lat[idx_max[2]]
# Ajuste de longitud para el mapa (-180 a 180)
lon_max_plot <- ifelse(lon_max > 180, lon_max - 360, lon_max)

# Encontrar coordenadas del Máximo Negativo (Cruz Amarilla)
min_val <- min(cor_map, na.rm=TRUE)
idx_min <- which(cor_map == min_val, arr.ind = TRUE)
lon_min <- lon[idx_min[1]]
lat_min <- lat[idx_min[2]]
# Ajuste de longitud para el mapa
lon_min_plot <- ifelse(lon_min > 180, lon_min - 360, lon_min)

# IMPRIMIR EN CONSOLA (Para que apuntes los datos)
cat("\n=========================================\n")
cat("      COORDENADAS DE LOS DIPOLOS\n")
cat("=========================================\n")
cat(sprintf("CRUZ VERDE (Máx Positivo): r = %.3f\n", max_val))
cat(sprintf(" -> Coordenadas: Lon %.1f, Lat %.1f\n", lon_max_plot, lat_max))
cat("    (Significado: Viento AQUÍ = Lluvia/Crecimiento)\n\n")

cat(sprintf("CRUZ AMARILLA (Máx Negativo): r = %.3f\n", min_val))
cat(sprintf(" -> Coordenadas: Lon %.1f, Lat %.1f\n", lon_min_plot, lat_min))
cat("    (Significado: Viento AQUÍ = Sequía/Decrecimiento)\n")
cat("=========================================\n")

# --- GRAFICAR ---
grid <- expand.grid(lon=lon, lat=lat)
grid$corr <- as.vector(cor_map)
grid$lon2 <- ifelse(grid$lon > 180, grid$lon - 360, grid$lon) # Ajuste mapa

world <- map_data("world")

ggplot() +
  geom_tile(data=grid, aes(x=lon2, y=lat, fill=corr)) +
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill=NA, color="black", size=0.2) +
  
  # AÑADIR LAS CRUCES
  geom_point(aes(x=lon_max_plot, y=lat_max), color="green", shape=3, size=5, stroke=2) +
  geom_point(aes(x=lon_min_plot, y=lat_min), color="yellow", shape=3, size=5, stroke=2) +
  
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, name="r (Pearson)") +
  coord_fixed(xlim=c(-100, 40), ylim=c(10, 80)) +
  labs(title = paste("Correlación de campo (invierno DJF) -", nivel_hpa, "hPa"),
       subtitle = "Verde: correl. pos. uwind-cronología | Amarillo: correl. neg. uwind-cronología",
       x="Longitud", y="Latitud") +
  theme_bw()
