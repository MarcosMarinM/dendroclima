# ==============================================================================
# 16_composite_anomalies_wind.R
# Author: Marcos Marín-Martín
# Date: 2026-02-11
# Description:
#   Composite Anomaly Maps (Extreme Growth Years).
#   Selects years of highest and lowest tree growth to generate composite maps
#   of wind anomalies, establishing the atmospheric conditions tailored to extremes.
# ==============================================================================

library(ncdf4); library(dplyr); library(ggplot2); library(maps)

# --- 1. CONFIGURACIÓN ---
file_nc   <- "PLACEHOLDER/path/to/uwnd.mon.mean.nc"
file_tree <- "PLACEHOLDER/path/to/chronology.txt"
nivel_hpa <- 300  # O el nivel que te haya dado mejor resultado (ej. 500)

# --- 2. CARGA DE DATOS ---
print("Cargando datos...")
nc <- nc_open(file_nc)
lev_idx <- which(abs(ncvar_get(nc, "level") - nivel_hpa) < 1)
# Leemos todo el mapa
uwnd <- ncvar_get(nc, "uwnd", start=c(1,1,lev_idx,1), count=c(-1,-1,1,-1))
lat <- ncvar_get(nc, "lat"); lon <- ncvar_get(nc, "lon")
dates <- as.Date(ncvar_get(nc, "time")/24, origin="1800-01-01")
nc_close(nc)

# Calcular Invierno (DJF) para cada punto del mapa
years <- as.numeric(format(dates, "%Y"))
months <- as.numeric(format(dates, "%m"))
unique_years <- unique(years)[-1]

# Matriz 3D [Lon, Lat, Año] con el promedio de invierno
print("Calculando inviernos...")
winter_wind <- array(NA, dim=c(length(lon), length(lat), length(unique_years)))

for(i in 1:length(unique_years)) {
  yr <- unique_years[i]
  idx <- which((years == (yr-1) & months == 12) | (years == yr & months %in% 1:2))
  if(length(idx) == 3) winter_wind[,,i] <- apply(uwnd[,,idx], c(1,2), mean)
}

# --- 3. SELECCIÓN DE AÑOS EXTREMOS ---
tree_data <- read.table(file_tree, header = TRUE)
# Filtramos solo años que tenemos en el clima
tree_sub <- tree_data %>% filter(year %in% unique_years)

# Umbral: Top 10% y Bottom 10% (o fijos, ej. 15 años)
n_extremes <- 11
high_years <- tree_sub %>% arrange(desc(res)) %>% head(n_extremes) %>% pull(year)
low_years  <- tree_sub %>% arrange(res) %>% head(n_extremes) %>% pull(year)

print(paste("Años de Crecimiento ALTO:", paste(high_years, collapse=", ")))
print(paste("Años de Crecimiento BAJO:", paste(low_years, collapse=", ")))

# --- 4. CÁLCULO DEL MAPA COMPUESTO ---
# Indices de los años
idx_high <- match(high_years, unique_years)
idx_low  <- match(low_years, unique_years)

# Promedios
map_high <- apply(winter_wind[,,idx_high], c(1,2), mean, na.rm=TRUE)
map_low  <- apply(winter_wind[,,idx_low], c(1,2), mean, na.rm=TRUE)

# DIFERENCIA (ALTO - BAJO)
# Esto nos dice: "¿Qué tiene de especial la atmósfera cuando el árbol crece mucho?"
map_diff <- map_high - map_low

# --- 5. GRAFICAR ---
grid <- expand.grid(lon=lon, lat=lat)
grid$diff <- as.vector(map_diff)
grid$lon2 <- ifelse(grid$lon > 180, grid$lon - 360, grid$lon)

world <- map_data("world")

ggplot() +
  geom_tile(data=grid, aes(x=lon2, y=lat, fill=diff)) +
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill=NA, color="black", size=0.2) +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, name="Wind anomaly (300 hPa) (m/s)") +
  coord_fixed(xlim=c(-100, 40), ylim=c(10, 80)) +
  labs(title = paste("Composite map: extreme growth years"),
       subtitle = paste("Difference in zonal wind (",nivel_hpa, "hPa). Red = stronger winds during high growth years"),
       x="Longitude", y="Latitude") +
  theme_bw()
