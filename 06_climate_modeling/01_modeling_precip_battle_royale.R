# ==============================================================================
# 01_modeling_precip_battle_royale.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Comparative Modeling of Precipitation using Atmospheric Variables ("Battle Royale").
#   Systematically compares linear models predicting precipitation using subsets 
#   of atmospheric variables: Zonal/Meridional Wind, Omega (Vertical Velocity), 
#   Surface Pressure, Relative Humidity, Precipitable Water, and Cloud Cover.
#
#   Methodology:
#   - Extracts climate data from NetCDF reanalysis files for a specific geographic box.
#   - Aggregates daily/monthly data into seasonal indices (e.g., Winter DJF).
#   - Fits multiple linear regression models (Univaraite, Bivariate, Trivariate).
#   - Evaluates models based on Adjusted R-squared (explanatory power) and AIC (efficiency).
#
#   Inputs:
#   - Monthly precipitation text file.
#   - NetCDF files for atmospheric variables (NCEP/NCAR Reanalysis).
#
#   Outputs:
#   - Console ranking of models.
#   - Lollipop chart visualizing model performance vs. complexity (R2 vs AIC).
# ==============================================================================

library(ncdf4)
library(dplyr)
library(ggplot2)
library(tidyr)

# ==============================================================================
# 1. PARAMETERS
# ==============================================================================
INPUT_PRECIP_FILE <- "PLACEHOLDER/path/to/precip_media.txt"
DIR_NETCDF        <- "PLACEHOLDER/path/to/20th_century_reanalysis/"

# Geometric definitions (Atlas Box)
LAT_BOX           <- c(30, 37)
LON_BOX           <- c(350, 5)

# NetCDF Filenames
FILE_UWND         <- file.path(DIR_NETCDF, "uwnd.mon.mean.nc")       # Zonal Wind
FILE_VWND         <- file.path(DIR_NETCDF, "vwnd.mon.mean.nc")       # Meridional Wind
FILE_OMEGA        <- file.path(DIR_NETCDF, "omega.mon.mean.nc")      # Vertical Velocity
FILE_PRMSL        <- file.path(DIR_NETCDF, "prmsl.mon.mean.nc")      # Surface Pressure
FILE_RHUM         <- file.path(DIR_NETCDF, "rhum.mon.mean.nc")       # Relative Humidity
FILE_PRWTR        <- file.path(DIR_NETCDF, "pr_wtr.eatm.mon.mean.nc")# Precipitable Water
FILE_TCDC         <- file.path(DIR_NETCDF, "tcdc.eatm.mon.mean.nc")  # Total Cloud Cover

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

#' Extract and aggregate level-based NetCDF data
get_level <- function(file, var, level, lat_b, lon_b) {
  nc <- nc_open(file)
  on.exit(nc_close(nc))
  
  levs <- ncvar_get(nc, "level")
  lev_idx <- which.min(abs(levs - level))
  
  lat <- ncvar_get(nc, "lat")
  lon <- ncvar_get(nc, "lon")
  time <- ncvar_get(nc, "time")
  dates <- as.Date(time/24, origin="1800-01-01")
  
  idx_lat <- which(lat >= lat_b[1] & lat <= lat_b[2])
  idx_lon <- c(which(lon >= lon_b[1]), which(lon <= lon_b[2]))
  
  if(length(idx_lat) == 0 || length(idx_lon) == 0) stop("No grid points found in specified box.")
  
  raw <- ncvar_get(nc, var, start=c(1,1,lev_idx,1), count=c(-1,-1,1,-1))
  val_box <- raw[idx_lon, idx_lat, ]
  
  # Average over space
  series <- apply(val_box, 3, mean, na.rm=TRUE)
  
  df <- data.frame(year=as.numeric(format(dates,"%Y")), month=as.numeric(format(dates,"%m")), val=series)
  df$clim_year <- ifelse(df$month==12, df$year+1, df$year)
  
  # DJF aggregation
  df %>% filter(month %in% c(12,1,2)) %>% 
    group_by(clim_year) %>% 
    summarise(mean=mean(val, na.rm=TRUE)) %>% 
    rename(year=clim_year)
}

#' Extract and aggregate surface NetCDF data
get_sfc <- function(file, var, lat_b, lon_b) {
  nc <- nc_open(file)
  on.exit(nc_close(nc))
  
  lat <- ncvar_get(nc, "lat")
  lon <- ncvar_get(nc, "lon")
  time <- ncvar_get(nc, "time")
  dates <- as.Date(time/24, origin="1800-01-01")
  
  idx_lat <- which(lat >= lat_b[1] & lat <= lat_b[2])
  idx_lon <- c(which(lon >= lon_b[1]), which(lon <= lon_b[2]))
  
  if(length(idx_lat) == 0 || length(idx_lon) == 0) stop("No grid points found in specified box.")
  
  raw <- ncvar_get(nc, var)
  val_box <- raw[idx_lon, idx_lat, ]
  
  # Average over space
  series <- apply(val_box, 3, mean, na.rm=TRUE)
  
  df <- data.frame(year=as.numeric(format(dates,"%Y")), month=as.numeric(format(dates,"%m")), val=series)
  df$clim_year <- ifelse(df$month==12, df$year+1, df$year)
  
  # DJF aggregation
  df %>% filter(month %in% c(12,1,2)) %>% 
    group_by(clim_year) %>% 
    summarise(mean=mean(val, na.rm=TRUE)) %>% 
    rename(year=clim_year)
}

# ==============================================================================
# 3. DATA LOADING
# ==============================================================================
message("--- CARGANDO VARIABLES FÍSICAS ---")

# Variable Y (Precipitation)
raw_p <- read.table(INPUT_PRECIP_FILE, header=TRUE)
raw_p$date_obj <- as.Date(as.character(raw_p$date), format="%Y%m%d")
raw_p$year <- as.numeric(format(raw_p$date_obj, "%Y"))
raw_p$month <- as.numeric(format(raw_p$date_obj, "%m"))
raw_p$clim_year <- ifelse(raw_p$month==12, raw_p$year+1, raw_p$year)

Y_dat <- raw_p %>% 
  filter(month %in% c(12,1,2)) %>% 
  group_by(clim_year) %>% 
  summarise(PCP=sum(precipitation, na.rm=TRUE), n=n()) %>% 
  filter(n>=3) %>% # Ensure complete season (approx)
  rename(year=clim_year)

# Variables X (Atmosphere)
message("1. Dinámica (Viento y Omega)...")
X_u300 <- get_level(FILE_UWND, "uwnd", 300, LAT_BOX, LON_BOX) %>% rename(U_300 = mean)
X_v300 <- get_level(FILE_VWND, "vwnd", 300, LAT_BOX, LON_BOX) %>% rename(V_300 = mean)
X_ome  <- get_level(FILE_OMEGA, "omega", 500, LAT_BOX, LON_BOX) %>% rename(Omega_500 = mean)

message("2. Termodinámica (humedad y nubes)...")
X_rh850 <- get_level(FILE_RHUM, "rhum", 850, LAT_BOX, LON_BOX) %>% rename(RH_850 = mean)
X_pw    <- get_sfc(FILE_PRWTR, "pr_wtr", LAT_BOX, LON_BOX) %>% rename(Pr_Wtr = mean)
X_cld   <- get_sfc(FILE_TCDC, "tcdc", LAT_BOX, LON_BOX) %>% rename(Clouds = mean)

message("3. Estado Sinóptico (Presión)...")
X_pres  <- get_sfc(FILE_PRMSL, "prmsl", LAT_BOX, LON_BOX) %>% rename(Press = mean)

# MERGE ALL DATA
dat <- Y_dat %>%
  inner_join(X_u300, by="year") %>% 
  inner_join(X_v300, by="year") %>%
  inner_join(X_ome, by="year")  %>% 
  inner_join(X_pres, by="year") %>%
  inner_join(X_rh850, by="year") %>% 
  inner_join(X_pw, by="year") %>%
  inner_join(X_cld, by="year")

# ==============================================================================
# 4. MODEL DEFINITIONS
# ==============================================================================
formulas <- list(
  "M01: Solo Viento (Base)"             = "PCP ~ U_300",
  "M02: Solo Presión"                   = "PCP ~ Press",
  "M03: Solo Humedad"                   = "PCP ~ RH_850",
  "M04: Solo Omega"                     = "PCP ~ Omega_500",
  "M05: Viento + Presión"               = "PCP ~ U_300 + Press",
  "M06: Viento + Humedad"               = "PCP ~ U_300 + RH_850",
  "M07: Viento + Omega"                 = "PCP ~ U_300 + Omega_500",
  "M08: Presión + Humedad"              = "PCP ~ Press + RH_850",
  "M09: PRO (Viento+Pres+Hum)"          = "PCP ~ U_300 + Press + RH_850",
  "M10: PRO + Nubes"                    = "PCP ~ U_300 + Press + RH_850 + Clouds",
  "M11: Flujo (Viento * Humedad)"       = "PCP ~ U_300 * RH_850",
  "M12: Vector Total (U + V + Pres)"    = "PCP ~ U_300 + V_300 + Press"
)

# ==============================================================================
# 5. MODEL COMPETITION (Ranking)
# ==============================================================================
results <- data.frame(Model = character(), R2_Adj = numeric(), AIC = numeric(), Num_Vars = numeric(), stringsAsFactors=FALSE)

for(name in names(formulas)) {
  mod <- lm(as.formula(formulas[[name]]), data = dat)
  s <- summary(mod)
  
  results[nrow(results)+1,] <- list(
    Model = name,
    R2_Adj = s$adj.r.squared * 100, # In percentage
    AIC = AIC(mod),
    Num_Vars = length(coef(mod)) - 1
  )
}

# Ranking: Prioritize Adjusted R2
ranking <- results %>% arrange(desc(R2_Adj))

# ==============================================================================
# 6. REPORTING
# ==============================================================================
cat("\n=======================================================\n")
cat("          🏆 RANKING DE MODELOS ATMOSFÉRICOS 🏆\n")
cat("=======================================================\n")
print(ranking)
cat("\n--- Interpretación ---\n")
cat("- R2_Adj: Cuanto más alto, mejor explica la lluvia (penalizando complejidad).\n")
cat("- AIC: Cuanto más BAJO, mejor equilibrio entre precisión y simplicidad.\n")
cat("=======================================================\n")

# ==============================================================================
# 7. VISUALIZATION (Leaderboard)
# ==============================================================================
results$Model <- factor(results$Model, levels = results$Model[order(results$R2_Adj)])

p <- ggplot(results, aes(x = R2_Adj, y = Model)) +
  geom_segment(aes(x = 0, xend = R2_Adj, y = Model, yend = Model), color = "grey") +
  geom_point(aes(color = as.factor(Num_Vars), size = AIC), alpha=0.9) +
  scale_color_brewer(palette = "Set1", name = "Nº Variables") +
  scale_size(range = c(5, 2), name = "AIC (Lower is better)") + 
  labs(title = "Comparativa de Modelos: ¿Qué explica la lluvia?",
       subtitle = "Eje X: Porcentaje de Lluvia Explicada (R2 Adj). Tamaño punto: Eficiencia.",
       x = "R2 Ajustado (%)", y = "") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, face = "bold"))

print(p)