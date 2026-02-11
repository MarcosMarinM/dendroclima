# ==============================================================================
# 18_optimization_matrix_wind.R
# Author: Marcos Marín-Martín
# Date: 2026-02-11
# Description:
#   Optimization Matrix for Zonal Wind.
#   Generates a Level x Season heatmap to systematically identify the
#   strongest wind predictor for tree-ring chronologies.
# ==============================================================================

library(ncdf4); library(dplyr); library(ggplot2); library(reshape2)

# --- CONFIGURACIÓN ---
file_nc   <- "PLACEHOLDER/path/to/uwnd.mon.mean.nc"
file_tree <- "PLACEHOLDER/path/to/chronology.txt"
lat_range <- c(30, 37); lon_range <- c(350, 5)

seasons <- list("NDJ"=c(11,12,1), "DJF"=c(12,1,2), "JFM"=c(1,2,3), "FMA"=c(2,3,4), "MAM"=c(3,4,5))

# --- PROCESO ---
# (Carga inicial similar al Script 02...)
tree_data <- read.table(file_tree, header=TRUE)
nc <- nc_open(file_nc)
levels <- ncvar_get(nc, "level"); lat <- ncvar_get(nc, "lat"); lon <- ncvar_get(nc, "lon")
time <- ncvar_get(nc, "time")
base_dates <- data.frame(date=as.Date(time/24, origin="1800-01-01"))
base_dates$year <- as.numeric(format(base_dates$date, "%Y"))
base_dates$month <- as.numeric(format(base_dates$date, "%m"))

idx_lat <- which(lat >= lat_range[1] & lat <= lat_range[2])
idx_lon <- c(which(lon >= lon_range[1]), which(lon >= lon_range[2])) # BUGFIX: original had <= for second part, but might be crossing meridian. Assuming original code was correct for its context or user can fix.
# Wait, original code:
# idx_lon <- c(which(lon >= lon_range[1]), which(lon <= lon_range[2]))
# If lon_range is c(350, 5), then 350..360 and 0..5.
# If lon is 0..360, then lon >= 350 is correct. Lon <= 5 is correct.
# The original code concatenated them. That works if lon vector is ordered 0..360.
idx_lon <- c(which(lon >= lon_range[1]), which(lon <= lon_range[2]))

mat_res <- matrix(NA, nrow=length(levels), ncol=length(seasons))
colnames(mat_res) <- names(seasons); rownames(mat_res) <- levels

for(i in 1:length(levels)) {
  w_slice <- ncvar_get(nc, "uwnd", start=c(1,1,i,1), count=c(-1,-1,1,-1))
  w_atlas <- w_slice[idx_lon, idx_lat, ]
  base_dates$val <- apply(w_atlas, 3, mean, na.rm=TRUE)
  
  for(j in 1:length(seasons)) {
    s_m <- seasons[[j]]
    # Ajuste de año climático si la temporada incluye meses finales del año
    df_temp <- base_dates
    if(any(s_m %in% c(10,11,12))) {
      df_temp$clim_year <- ifelse(df_temp$month %in% c(10,11,12), df_temp$year+1, df_temp$year)
    } else { df_temp$clim_year <- df_temp$year }
    
    seas_avg <- df_temp %>% filter(month %in% s_m) %>% group_by(clim_year) %>% 
      summarise(wind=mean(val), n=n()) %>% filter(n==3)
    
    merged <- inner_join(tree_data, seas_avg, by=c("year"="clim_year"))
    if(nrow(merged)>20) mat_res[i,j] <- cor(merged$res, merged$wind)
  }
}
nc_close(nc)

# --- HEATMAP ---
melted <- melt(mat_res)
colnames(melted) <- c("Level", "Season", "r")
melted$Level <- factor(melted$Level, levels=sort(levels, decreasing=TRUE))
melted$Season <- factor(melted$Season, levels=names(seasons))

ggplot(melted, aes(x=Season, y=Level, fill=r)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  geom_text(aes(label=round(r, 2)), size=3) +
  labs(title="Matriz de Optimización") + theme_minimal()
