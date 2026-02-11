# ==============================================================================
# 19_vertical_profile_correlation.R
# Author: Marcos Marín-Martín
# Date: 2026-02-11
# Description:
#   Vertical Profile of Correlations.
#   Plots the correlation coefficient (X-axis) against Geopotential Height/Pressure (Y-axis)
#   to visualize the vertical structure of the climate signal in tree rings.
# ==============================================================================

library(ncdf4); library(dplyr); library(ggplot2)

# --- CONFIGURACIÓN ---
file_nc   <- "PLACEHOLDER/path/to/uwnd.mon.mean.nc"
file_tree <- "PLACEHOLDER/path/to/chronology.txt"

# Coordenadas Caja Atlas
lat_range <- c(30, 37); lon_range <- c(350, 5)

# --- PROCESO ---
tree_data <- read.table(file_tree, header = TRUE)
nc <- nc_open(file_nc)
levels <- ncvar_get(nc, "level"); lat <- ncvar_get(nc, "lat"); lon <- ncvar_get(nc, "lon")
time <- ncvar_get(nc, "time")
dates <- data.frame(date = as.Date(time/24, origin="1800-01-01"))
dates$year <- as.numeric(format(dates$date, "%Y"))
dates$month <- as.numeric(format(dates$date, "%m"))
dates$clim_year <- ifelse(dates$month == 12, dates$year + 1, dates$year)

# Indices espaciales
idx_lat <- which(lat >= lat_range[1] & lat <= lat_range[2])
idx_lon <- c(which(lon >= lon_range[1]), which(lon <= lon_range[2]))

results <- data.frame(level=numeric(), r=numeric(), p=numeric())

for(l in levels) {
  # Leer nivel l
  lev_idx <- which(abs(levels - l) < 1)
  w_slice <- ncvar_get(nc, "uwnd", start=c(1,1,lev_idx,1), count=c(-1,-1,1,-1))
  
  # Recortar y promediar espacialmente
  w_atlas <- w_slice[idx_lon, idx_lat, ]
  dates$val <- apply(w_atlas, 3, mean, na.rm=TRUE)
  
  # Calcular DJF
  djf_mean <- dates %>% filter(month %in% c(12,1,2)) %>% 
    group_by(clim_year) %>% summarise(wind=mean(val))
  
  merged <- inner_join(tree_data, djf_mean, by=c("year"="clim_year"))
  
  if(nrow(merged)>15) {
    test <- cor.test(merged$res, merged$wind)
    results <- rbind(results, data.frame(level=l, r=test$estimate, p=test$p.value))
  }
}
nc_close(nc)

# --- GRAFICAR ---
ggplot(results, aes(x=r, y=level)) +
  geom_path(color="gray") + geom_point(aes(color=p<0.05), size=3) +
  scale_y_reverse(breaks=levels) + # Invertir eje Y (atmosférico)
  geom_vline(xintercept=0, linetype="dashed") +
  labs(title="Perfil Vertical de Correlación (Atlas Box - DJF)", x="Pearson r", y="Presión (hPa)") +
  theme_minimal()
