# ==============================================================================
# 03_validation_model.R
# Author: Marcos Marín-Martín
# Date: 2026-02-12
# Description:
#   Generic Validation & Reconstruction Script.
#   - Takes ANY linear model formula as input.
#   - Loads necessary data.
#   - Exports statistical diagnostics to a specific output folder.
#   - Performs Leave-One-Out Cross-Validation (LOOCV).
# ==============================================================================

# LIBRARIES
library(ncdf4)
library(ggplot2)
library(tidyr)
library(car)    # VIF
library(lmtest) # Durbin-Watson
library(dplyr)  # Loaded last

# ==============================================================================
# 1. PARAMETERS & INPUTS
# ==============================================================================

# --- [IMPORTANTE] PEGA AQUÍ LA FÓRMULA GANADORA DEL STEP 02 ---
WINNING_FORMULA_STR <- "PCP ~ xxx + xxx + xxx"  # Ejemplo: "PCP ~ MOI + NAO + AMO + TSI + AOD + U_Atlas + Press + RH_850 + Omega + Clouds"

# --- CONFIGURACIÓN DE SALIDA (NUEVO) ---
# Define dónde quieres que se guarden el .txt y el .csv
OUTPUT_DIR <- "./results/validation/" 
OUTPUT_TXT_NAME <- "03_model_diagnostics.txt"
OUTPUT_CSV_NAME <- "03_reconstruction_data.csv"

# --- RUTAS DE ARCHIVOS (PLACEHOLDERS) ---
INPUT_PRECIP_FILE <- "PLACEHOLDER/path/to/precip_media.txt"
INPUT_MOI_FILE    <- "PLACEHOLDER/path/to/moi_index.txt"
INPUT_NAO_FILE    <- "PLACEHOLDER/path/to/nao_index.txt"
INPUT_AMO_FILE    <- "PLACEHOLDER/path/to/amo_index.txt"
INPUT_SOLAR_FILE  <- "PLACEHOLDER/path/to/solar_forcing.txt"
INPUT_VOLC_FILE   <- "PLACEHOLDER/path/to/volc_forcing.txt"

# --- REANALYSIS DATA ---
DIR_NETCDF        <- "PLACEHOLDER/path/to/20th_century_reanalysis/"

# --- GEOMETRIC DEFINITIONS ---
BOX_ATLAS   <- list(lat=c(30, 37), lon=c(350, 5))
BOX_IRELAND <- list(lat=c(50, 60), lon=c(340, 360))

# --- NETCDF FILENAMES ---
FILE_UWND   <- file.path(DIR_NETCDF, "uwnd.mon.mean.nc")
FILE_VWND   <- file.path(DIR_NETCDF, "vwnd.mon.mean.nc")
FILE_OMEGA  <- file.path(DIR_NETCDF, "omega.mon.mean.nc")
FILE_PRMSL  <- file.path(DIR_NETCDF, "prmsl.mon.mean.nc")
FILE_RHUM   <- file.path(DIR_NETCDF, "rhum.mon.mean.nc")
FILE_PRWTR  <- file.path(DIR_NETCDF, "pr_wtr.eatm.mon.mean.nc") 
FILE_TCDC   <- file.path(DIR_NETCDF, "tcdc.eatm.mon.mean.nc")

# Crear directorio de salida si no existe
if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive=TRUE)

# ==============================================================================
# 2. ROBUST DATA LOADING
# ==============================================================================
get_level_val <- function(file, var, level, box) {
  if(!file.exists(file)) return(NULL)
  nc <- nc_open(file); on.exit(nc_close(nc))
  levs <- ncvar_get(nc, "level")
  if(!is.null(levs)) {
    lev_idx <- which.min(abs(levs - level))
    start_v <- c(1,1,lev_idx,1); count_v <- c(-1,-1,1,-1)
  } else { start_v <- c(1,1,1); count_v <- c(-1,-1,-1) }
  
  lat <- ncvar_get(nc,"lat"); lon <- ncvar_get(nc,"lon")
  time <- ncvar_get(nc,"time"); dates <- as.Date(time/24, origin="1800-01-01")
  
  idx_lat <- which(lat >= box$lat[1] & lat <= box$lat[2])
  idx_lon <- c(which(lon >= box$lon[1]), which(lon <= box$lon[2]))
  
  raw <- tryCatch(ncvar_get(nc, var, start=start_v, count=count_v), error=function(e) NULL)
  if(is.null(raw)) return(NULL)
  
  val_box <- raw[idx_lon, idx_lat, ]
  series <- apply(val_box, 3, mean, na.rm=TRUE)
  
  df <- data.frame(year=as.numeric(format(dates,"%Y")), month=as.numeric(format(dates,"%m")), val=series)
  df$clim_year <- ifelse(df$month==12, df$year+1, df$year)
  df %>% filter(month %in% c(12,1,2)) %>% group_by(clim_year) %>% summarise(mean=mean(val, na.rm=TRUE)) %>% rename(year=clim_year)
}

message("\n--- 1. CARGANDO DATOS ---")
data_list <- list()

# 2.1 PRECIPITATION
if(file.exists(INPUT_PRECIP_FILE)) {
  raw_p <- read.table(INPUT_PRECIP_FILE, header=TRUE)
  raw_p$date_obj <- as.Date(as.character(raw_p$date), format="%Y%m%d")
  raw_p$year <- as.numeric(format(raw_p$date_obj, "%Y"))
  raw_p$month <- as.numeric(format(raw_p$date_obj, "%m"))
  raw_p$clim_year <- ifelse(raw_p$month==12, raw_p$year+1, raw_p$year)
  data_list[["PCP"]] <- raw_p %>% filter(month %in% c(12,1,2)) %>% 
    group_by(clim_year) %>% summarise(PCP=sum(precipitation, na.rm=T), n=n()) %>% 
    filter(n>=3) %>% select(clim_year, PCP) %>% rename(year=clim_year)
} else { stop("❌ ERROR: Precip file missing.") }

# 2.2 INDICES
load_idx <- function(path, name, col_name) {
  if(file.exists(path)) {
    df <- read.table(path, header=TRUE)
    if("date" %in% names(df)) { # MOI case
      df$date_obj <- as.Date(as.character(df$date), format="%Y%m%d")
      df$year <- as.numeric(format(df$date_obj, "%Y")); df$month <- as.numeric(format(df$date_obj, "%m"))
      df$clim_year <- ifelse(df$month==12, df$year+1, df$year)
      df <- df %>% filter(month %in% c(12,1,2)) %>% group_by(clim_year) %>% 
        summarise(val=mean(get(col_name), na.rm=T)) %>% rename(year=clim_year)
    } else { 
      if("Year" %in% names(df)) df <- df %>% rename(year=Year)
      df <- df %>% select(year, all_of(col_name)) %>% rename(val=all_of(col_name))
    }
    names(df)[2] <- name
    return(df)
  }
  return(NULL)
}

data_list[["MOI"]]   <- load_idx(INPUT_MOI_FILE, "MOI", "moi")
data_list[["NAO"]]   <- load_idx(INPUT_NAO_FILE, "NAO", "NAOmed")
data_list[["AMO"]]   <- load_idx(INPUT_AMO_FILE, "AMO", "AMO")
data_list[["TSI"]]   <- load_idx(INPUT_SOLAR_FILE, "TSI", "TSI")
data_list[["AOD"]]   <- load_idx(INPUT_VOLC_FILE, "AOD", "AOD")

# 2.3 REANALYSIS VARS
if(dir.exists(DIR_NETCDF)) {
  tmp <- get_level_val(FILE_UWND, "uwnd", 300, BOX_ATLAS); if(!is.null(tmp)) data_list[["U_Atlas"]] <- tmp %>% rename(U_Atlas=mean)
  tmp <- get_level_val(FILE_PRMSL, "prmsl", -999, BOX_ATLAS); if(!is.null(tmp)) data_list[["Press"]] <- tmp %>% rename(Press=mean)
  tmp <- get_level_val(FILE_RHUM, "rhum", 850, BOX_ATLAS); if(!is.null(tmp)) data_list[["RH_850"]] <- tmp %>% rename(RH_850=mean)
  tmp <- get_level_val(FILE_OMEGA, "omega", 500, BOX_ATLAS); if(!is.null(tmp)) data_list[["Omega"]] <- tmp %>% rename(Omega=mean)
  tmp <- get_level_val(FILE_TCDC, "tcdc", -999, BOX_ATLAS); if(!is.null(tmp)) data_list[["Clouds"]] <- tmp %>% rename(Clouds=mean)
}

dat <- Reduce(function(x, y) inner_join(x, y, by="year"), data_list)
dat <- na.omit(dat)
message(paste("✅ Datos consolidados:", nrow(dat), "años."))

# ==============================================================================
# 3. MODEL CALIBRATION & TXT EXPORT
# ==============================================================================
message("\n--- 2. CALIBRANDO Y EXPORTANDO DIAGNÓSTICOS ---")

final_formula <- as.formula(WINNING_FORMULA_STR)
model <- lm(final_formula, data=dat)

# Construir ruta completa del archivo de salida
out_txt_path <- file.path(OUTPUT_DIR, OUTPUT_TXT_NAME)

sink(out_txt_path) # Redirigir salida al archivo en la carpeta definida

cat("==============================================================\n")
cat(" CLIMATE RECONSTRUCTION MODEL DIAGNOSTICS\n")
cat(sprintf(" Date: %s\n", Sys.Date()))
cat(sprintf(" Formula: %s\n", WINNING_FORMULA_STR))
cat("==============================================================\n\n")

cat("1. MODEL SUMMARY\n")
cat("----------------\n")
print(summary(model))

cat("\n2. MULTICOLLINEARITY (VIF)\n")
cat("--------------------------\n")
vif_val <- tryCatch(vif(model), error=function(e) "Not applicable (Univariate)")
print(vif_val)

cat("\n3. AUTOCORRELATION (Durbin-Watson)\n")
cat("----------------------------------\n")
print(dwtest(model))

cat("\n4. RESIDUAL NORMALITY (Shapiro-Wilk)\n")
cat("------------------------------------\n")
print(shapiro.test(resid(model)))

cat("\n5. COEFFICIENTS (For Equation)\n")
cat("------------------------------\n")
print(coef(model))

sink() # Stop writing
message(paste("✅ Diagnósticos guardados en:", out_txt_path))

# ==============================================================================
# 4. CROSS-VALIDATION (LOOCV) GENERIC
# ==============================================================================
message("\n--- 3. EJECUTANDO VALIDACIÓN CRUZADA (LOOCV) ---")

cv_preds <- numeric(nrow(dat))

for(i in 1:nrow(dat)) {
  train <- dat[-i, ]
  test  <- dat[i, ]
  m_tmp <- lm(final_formula, data=train)
  cv_preds[i] <- predict(m_tmp, newdata=test)
}

dat$PCP_CV <- cv_preds
sse <- sum((dat$PCP - dat$PCP_CV)^2)
sst <- sum((dat$PCP - mean(dat$PCP))^2)
r2_val <- 1 - (sse/sst)
r2_cal <- summary(model)$adj.r.squared

cat(sprintf(">>> R2 Calibración: %.2f%%\n", r2_cal*100))
cat(sprintf(">>> R2 Validación:  %.2f%%\n", r2_val*100))

# ==============================================================================
# 5. PLOTS & EXPORT
# ==============================================================================

# Gráfico Serie Temporal
p1 <- ggplot(dat, aes(x=year)) +
  geom_line(aes(y=PCP, color="Observed"), linewidth=0.8) +
  geom_line(aes(y=PCP_CV, color="Reconstructed (CV)"), linewidth=0.8, linetype="solid") +
  scale_color_manual(values=c("Observed"="black", "Reconstructed (CV)"="red")) +
  labs(title="Reconstruction vs Observation",
       subtitle=paste0("Formula: ", WINNING_FORMULA_STR, "\nR2 Val: ", round(r2_val*100,2), "%"),
       y="Precipitation (mm)", x="Year", color="") +
  theme_minimal() + theme(legend.position="bottom")

# Gráfico Scatter
p2 <- ggplot(dat, aes(x=PCP, y=PCP_CV)) +
  geom_point(alpha=0.6) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey") +
  geom_smooth(method="lm", se=FALSE, color="red", linewidth=0.5) +
  labs(title="Goodness of Fit", x="Observed", y="Predicted (LOOCV)") +
  theme_minimal()

print(p1)
print(p2)

# Guardar CSV en la carpeta definida
out_csv_path <- file.path(OUTPUT_DIR, OUTPUT_CSV_NAME)
write.csv(dat, out_csv_path, row.names=FALSE)
message(paste("✅ Datos reconstruidos guardados en:", out_csv_path))