# ==============================================================================
# 02_modeling_precip_ultimate_battle.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   Comprehensive Climate Modeling "Ultimate Battle" (19 Models).
#   An extensive comparison of 19 statistical models for predicting precipitation,
#   integrating local atmospheric dynamics (Wind, Pressure, Humidity) with 
#   large-scale teleconnection indices (NAO, AMO, Western Mediterranean Oscillation Dipole).
#
#   Methodology:
#   - Integrates separate data sources: Precipitation (Station), Indices (NAO, AMO), and Reanalysis (NetCDF).
#   - Calculates a custom dipole index (Atlas vs Ireland Zonal Wind).
#   - Tests nested models ranging from simple univariate regressions to complex multivariable combinations.
#   - Generates explicit regression equations for each model.
#
#   Inputs:
#   - Precipitation data.
#   - Teleconnection indices (NAO, AMO).
#   - NetCDF atmospheric data.
#
#   Outputs:
#   - Detailed ranking table (R2 adj, AIC).
#   - Specific regression equations.
#   - "Lollipop" efficiency plot.
# ==============================================================================

library(ncdf4)
library(dplyr)
library(ggplot2)
library(tidyr)
library(zoo)

# ==============================================================================
# 1. PARAMETERS
# ==============================================================================
INPUT_PRECIP_FILE <- "PLACEHOLDER/path/to/precip_media.txt"
INPUT_NAO_FILE    <- "PLACEHOLDER/path/to/nao_index.txt"
INPUT_AMO_FILE    <- "PLACEHOLDER/path/to/amo_index.txt"
DIR_NETCDF        <- "PLACEHOLDER/path/to/20th_century_reanalysis/"

# Box Definitions
BOX_ATLAS   <- list(lat=c(30, 37), lon=c(350, 5))
BOX_IRELAND <- list(lat=c(50, 60), lon=c(340, 360))

# NetCDF Filenames
FILE_UWND   <- file.path(DIR_NETCDF, "uwnd.mon.mean.nc")
FILE_VWND   <- file.path(DIR_NETCDF, "vwnd.mon.mean.nc")
FILE_OMEGA  <- file.path(DIR_NETCDF, "omega.mon.mean.nc")
FILE_PRMSL  <- file.path(DIR_NETCDF, "prmsl.mon.mean.nc")
FILE_RHUM   <- file.path(DIR_NETCDF, "rhum.mon.mean.nc")        # Levels (850hPa)
FILE_PRWTR  <- file.path(DIR_NETCDF, "pr_wtr.eatm.mon.mean.nc") 
FILE_TCDC   <- file.path(DIR_NETCDF, "tcdc.eatm.mon.mean.nc")

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

get_level <- function(file, var, level, box) {
  nc <- nc_open(file)
  on.exit(nc_close(nc))
  
  levs <- ncvar_get(nc, "level"); lev_idx <- which.min(abs(levs - level))
  
  lat <- ncvar_get(nc,"lat"); lon <- ncvar_get(nc,"lon")
  time <- ncvar_get(nc,"time"); dates <- as.Date(time/24, origin="1800-01-01")
  
  idx_lat <- which(lat >= box$lat[1] & lat <= box$lat[2])
  idx_lon <- c(which(lon >= box$lon[1]), which(lon <= box$lon[2]))
  
  raw <- ncvar_get(nc, var, start=c(1,1,lev_idx,1), count=c(-1,-1,1,-1))
  val_box <- raw[idx_lon, idx_lat, ]
  
  series <- apply(val_box, 3, mean, na.rm=TRUE)
  
  df <- data.frame(year=as.numeric(format(dates,"%Y")), month=as.numeric(format(dates,"%m")), val=series)
  df$clim_year <- ifelse(df$month==12, df$year+1, df$year)
  
  df %>% filter(month %in% c(12,1,2)) %>% 
    group_by(clim_year) %>% 
    summarise(mean=mean(val, na.rm=TRUE)) %>% 
    rename(year=clim_year)
}

get_sfc <- function(file, var, box) {
  nc <- nc_open(file)
  on.exit(nc_close(nc))
  
  lat <- ncvar_get(nc,"lat"); lon <- ncvar_get(nc,"lon")
  time <- ncvar_get(nc,"time"); dates <- as.Date(time/24, origin="1800-01-01")
  
  idx_lat <- which(lat >= box$lat[1] & lat <= box$lat[2])
  idx_lon <- c(which(lon >= box$lon[1]), which(lon <= box$lon[2]))
  
  raw <- ncvar_get(nc, var)
  val_box <- raw[idx_lon, idx_lat, ]
  
  series <- apply(val_box, 3, mean, na.rm=TRUE)
  
  df <- data.frame(year=as.numeric(format(dates,"%Y")), month=as.numeric(format(dates,"%m")), val=series)
  df$clim_year <- ifelse(df$month==12, df$year+1, df$year)
  
  df %>% filter(month %in% c(12,1,2)) %>% 
    group_by(clim_year) %>% 
    summarise(mean=mean(val, na.rm=TRUE)) %>% 
    rename(year=clim_year)
}

# ==============================================================================
# 3. DATA LOADING
# ==============================================================================
message("--- CARGANDO VARIABLES... ---")

# Precipitation
raw_p <- read.table(INPUT_PRECIP_FILE, header=TRUE)
raw_p$date_obj <- as.Date(as.character(raw_p$date), format="%Y%m%d")
raw_p$year <- as.numeric(format(raw_p$date_obj, "%Y"))
raw_p$month <- as.numeric(format(raw_p$date_obj, "%m"))
raw_p$clim_year <- ifelse(raw_p$month==12, raw_p$year+1, raw_p$year)

Y_dat <- raw_p %>% 
  filter(month %in% c(12,1,2)) %>% 
  group_by(clim_year) %>% 
  summarise(PCP=sum(precipitation, na.rm=TRUE), n=n()) %>% 
  filter(n>=3) %>% 
  rename(year=clim_year)

# Local Atmospheric Variables
X_u300  <- get_level(FILE_UWND, "uwnd", 300, BOX_ATLAS) %>% rename(U_Atlas = mean)
X_v300  <- get_level(FILE_VWND, "vwnd", 300, BOX_ATLAS) %>% rename(V_Atlas = mean)
X_ome   <- get_level(FILE_OMEGA, "omega", 500, BOX_ATLAS) %>% rename(Omega = mean)
X_rh850 <- get_level(FILE_RHUM, "rhum", 850, BOX_ATLAS) %>% rename(RH_850 = mean)
X_pres  <- get_sfc(FILE_PRMSL, "prmsl", BOX_ATLAS) %>% rename(Press = mean)
X_prw   <- get_sfc(FILE_PRWTR, "pr_wtr", BOX_ATLAS) %>% rename(Pr_Wtr = mean)
X_cld   <- get_sfc(FILE_TCDC, "tcdc", BOX_ATLAS) %>% rename(Clouds = mean)
X_uIre  <- get_level(FILE_UWND, "uwnd", 300, BOX_IRELAND) %>% rename(U_Ireland = mean)

# Calculate Dipole Index
X_dip   <- inner_join(X_u300, X_uIre, by="year") %>% 
  mutate(DIP_Index = U_Atlas - U_Ireland) %>% 
  select(year, DIP_Index)

# NAO
nao_raw <- read.table(INPUT_NAO_FILE, header=TRUE)
X_nao   <- nao_raw %>% select(Year, NAOmed) %>% rename(year = Year, NAO = NAOmed)

# AMO (Winter Calculation)
amo_raw <- read.table(INPUT_AMO_FILE, header=FALSE)
colnames(amo_raw) <- c("Year", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
amo_long <- amo_raw %>%
  pivot_longer(cols = -Year, names_to = "Month_Name", values_to = "Value") %>%
  mutate(Month = match(Month_Name, month.abb)) %>%
  mutate(clim_year = ifelse(Month == 12, Year + 1, Year)) %>%
  filter(Month %in% c(12, 1, 2)) %>%
  group_by(clim_year) %>%
  summarise(AMO = mean(Value, na.rm=TRUE), n=n()) %>%
  filter(n>=3) %>% 
  rename(year = clim_year)

# MASTER JOIN
dat <- Y_dat %>%
  inner_join(X_u300, by="year") %>% 
  inner_join(X_v300, by="year") %>%
  inner_join(X_ome, by="year")  %>% 
  inner_join(X_pres, by="year") %>%
  inner_join(X_rh850, by="year")%>% 
  inner_join(X_prw, by="year")  %>%
  inner_join(X_cld, by="year")  %>% 
  inner_join(X_dip, by="year")  %>%
  inner_join(X_nao, by="year")  %>% 
  inner_join(amo_long, by="year")

# ==============================================================================
# 4. MODEL LIST
# ==============================================================================
formulas <- list(
  # --- UNIVARIATE (SIMPLE) ---
  "M01_Wind"        = "PCP ~ U_Atlas",
  "M02_Pres"        = "PCP ~ Press",
  "M03_Hum850"      = "PCP ~ RH_850",
  "M04_Omega"       = "PCP ~ Omega",
  "M05_PrWtr"       = "PCP ~ Pr_Wtr",
  "M06_Clouds"      = "PCP ~ Clouds",
  "M07_NAO"         = "PCP ~ NAO",
  "M08_AMO"         = "PCP ~ AMO",
  "M09_Dipolo"      = "PCP ~ DIP_Index",
  
  # --- BIVARIATE (MIXED) ---
  "M10_Wind_Pres"   = "PCP ~ U_Atlas + Press",
  "M11_Wind_Hum"    = "PCP ~ U_Atlas + RH_850",
  "M12_Pres_Hum"    = "PCP ~ Press + RH_850",
  "M13_Sup_Base"    = "PCP ~ NAO + DIP_Index",
  
  # --- TRIVARIATE (PHYSICS) ---
  "M14_PRO_Local"   = "PCP ~ U_Atlas + Press + RH_850",
  "M15_Vector3D"    = "PCP ~ U_Atlas + V_Atlas + Press",
  
  # --- MULTIVARIATE (COMBOS) ---
  "M16_PRO_Clouds"  = "PCP ~ U_Atlas + Press + RH_850 + Clouds",
  "M17_Combo_NAO"   = "PCP ~ NAO + U_Atlas + Press + RH_850",
  "M18_Combo_AMO"   = "PCP ~ U_Atlas + Press + RH_850 + AMO",
  "M19_ALL_STARS"   = "PCP ~ NAO + AMO + U_Atlas + Press + RH_850"
)

# ==============================================================================
# 5. EXECUTION LOOP
# ==============================================================================
ranking_df <- data.frame(Model=character(), R2_Adj=numeric(), AIC=numeric(), Num_Vars=numeric(), stringsAsFactors=FALSE)

cat("\n\n################################################################\n")
cat("   📚 ENCICLOPEDIA DE MODELOS (M01 - M19)\n")
cat("################################################################\n\n")

for(name in names(formulas)) {
  mod <- lm(as.formula(formulas[[name]]), data = dat)
  s <- summary(mod)
  cf <- coef(mod)
  
  ranking_df[nrow(ranking_df)+1,] <- list(name, s$adj.r.squared*100, AIC(mod), length(cf)-1)
  
  cat(paste0(">>> ", name, "\n"))
  cat(sprintf("   R2 Adj: %.2f%% | AIC: %.1f\n", s$adj.r.squared*100, AIC(mod)))
  
  # Generate Equation String
  eq_str <- sprintf("PCP = %.2f", cf[1])
  for(i in 2:length(cf)) {
    signo <- ifelse(cf[i] >= 0, "+", "-")
    val   <- abs(cf[i])
    var   <- names(cf)[i]
    eq_str <- paste0(eq_str, sprintf(" %s %.4f·%s", signo, val, var))
  }
  cat(paste0("   Eq: ", eq_str, "\n"))
  cat("----------------------------------------------------------------\n")
}

# ==============================================================================
# 6. VISUALIZATION
# ==============================================================================
ranking_df$Model <- factor(ranking_df$Model, levels=ranking_df$Model[order(ranking_df$R2_Adj)])

# Scale AIC for size (Inverted: Lower AIC = Larger Dot)
min_aic <- min(ranking_df$AIC); max_aic <- max(ranking_df$AIC)
ranking_df$Efficiency <- (max_aic - ranking_df$AIC) + 10 # Visual offset

p <- ggplot(ranking_df, aes(x=R2_Adj, y=Model)) +
  geom_segment(aes(x=0, xend=R2_Adj, y=Model, yend=Model), color="grey60", linetype="dotted") +
  geom_point(aes(color=as.factor(Num_Vars), size=Efficiency), alpha=0.9) +
  
  geom_text(aes(label=sprintf("%.1f%%", R2_Adj)), hjust=-0.3, size=3.5, fontface="bold") +
  
  scale_color_viridis_d(option = "plasma", name="N.º variables", end=0.9) +
  scale_size(range = c(2, 8), guide="none") + 
  
  labs(title="Batalla de modelos",
       subtitle="Más a la derecha = Mejor predicción. Bola más grande = modelo más eficiente (AIC).",
       x="R2 ajustado (%)", y="") +
  
  theme_minimal() +
  theme(axis.text.y = element_text(face="bold", size=9),
        panel.grid.major.y = element_blank()) +
  xlim(0, max(ranking_df$R2_Adj) + 8)

# Print Final Ranking
print(ranking_df %>% select(Model, R2_Adj, AIC, Num_Vars) %>% arrange(desc(R2_Adj)))
print(p)