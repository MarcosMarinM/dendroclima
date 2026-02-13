# ==============================================================================
# 02_modeling_precip_ultimate_battle.R
# Author: Marcos Marín-Martín
# Date: 2026-02-12
# Description:
#   Advanced climate modeling: linear, interactions, stepwise & GAMs.
#   Robust to missing data. 
#   Integrates: 
#     - Station data (precip)
#     - Teleconnections (NAO, AMO, MOI)
#     - Reanalysis (NetCDF)
#     - External forcing (solar, volcanoes)
# ==============================================================================

library(ncdf4)
library(ggplot2)
library(tidyr)
library(zoo)
library(mgcv)
library(MASS)   
library(dplyr)

# ==============================================================================
# 1. PARAMETERS & PATHS
# ==============================================================================
# --- INPUT FILES (STATION & INDICES) ---
INPUT_PRECIP_FILE <- "PLACEHOLDER/path/to/precip_media.txt"
INPUT_MOI_FILE    <- "PLACEHOLDER/path/to/moi_index.txt"
INPUT_NAO_FILE    <- "PLACEHOLDER/path/to/nao_index.txt"
INPUT_AMO_FILE    <- "PLACEHOLDER/path/to/amo_index.txt"
INPUT_SOLAR_FILE  <- "PLACEHOLDER/path/to/solar_forcing.txt"
INPUT_VOLC_FILE   <- "PLACEHOLDER/path/to/volc_forcing.txt"

# --- REANALYSIS DATA (NETCDF DIRECTORY) ---
DIR_NETCDF        <- "PLACEHOLDER/path/to/20th_century_reanalysis/"

# --- GEOMETRIC DEFINITIONS ---
BOX_ATLAS   <- list(lat=c(30, 37), lon=c(350, 5))
BOX_IRELAND <- list(lat=c(50, 60), lon=c(340, 360))

# --- NETCDF FILENAMES ---
FILE_UWND   <- file.path(DIR_NETCDF, "uwnd.mon.mean.nc")        # Zonal Wind
FILE_VWND   <- file.path(DIR_NETCDF, "vwnd.mon.mean.nc")        # Meridional Wind
FILE_OMEGA  <- file.path(DIR_NETCDF, "omega.mon.mean.nc")       # Vertical Velocity
FILE_PRMSL  <- file.path(DIR_NETCDF, "prmsl.mon.mean.nc")       # Surface Pressure
FILE_RHUM   <- file.path(DIR_NETCDF, "rhum.mon.mean.nc")        # Relative Humidity
FILE_PRWTR  <- file.path(DIR_NETCDF, "pr_wtr.eatm.mon.mean.nc") # Precipitable Water
FILE_TCDC   <- file.path(DIR_NETCDF, "tcdc.eatm.mon.mean.nc")   # Total Cloud Cover

# ==============================================================================
# 2. HELPER FUNCTIONS (NETCDF EXTRACTION)
# ==============================================================================

get_level <- function(file, var, level, box) {
  if(!file.exists(file)) return(NULL)
  nc <- nc_open(file)
  on.exit(nc_close(nc))
  
  levs <- ncvar_get(nc, "level")
  # Si hay niveles, buscar el más cercano. Si no, ignorar.
  if(!is.null(levs)) {
    lev_idx <- which.min(abs(levs - level))
    start_vec <- c(1,1,lev_idx,1)
    count_vec <- c(-1,-1,1,-1)
  } else {
    start_vec <- c(1,1,1)
    count_vec <- c(-1,-1,-1)
  }
  
  lat <- ncvar_get(nc,"lat"); lon <- ncvar_get(nc,"lon")
  time <- ncvar_get(nc,"time"); dates <- as.Date(time/24, origin="1800-01-01")
  
  idx_lat <- which(lat >= box$lat[1] & lat <= box$lat[2])
  idx_lon <- c(which(lon >= box$lon[1]), which(lon <= box$lon[2]))
  
  # Carga segura dependiendo de dimensiones
  raw <- tryCatch({
    ncvar_get(nc, var, start=start_vec, count=count_vec)
  }, error = function(e) return(NULL))
  
  if(is.null(raw)) return(NULL)
  
  val_box <- raw[idx_lon, idx_lat, ]
  series <- apply(val_box, 3, mean, na.rm=TRUE)
  
  df <- data.frame(year=as.numeric(format(dates,"%Y")), month=as.numeric(format(dates,"%m")), val=series)
  df$clim_year <- ifelse(df$month==12, df$year+1, df$year)
  
  df %>% filter(month %in% c(12,1,2)) %>% 
    group_by(clim_year) %>% summarise(mean=mean(val, na.rm=TRUE)) %>% rename(year=clim_year)
}

get_sfc <- function(file, var, box) {
  if(!file.exists(file)) return(NULL)
  nc <- nc_open(file)
  on.exit(nc_close(nc))
  
  lat <- ncvar_get(nc,"lat"); lon <- ncvar_get(nc,"lon")
  time <- ncvar_get(nc,"time"); dates <- as.Date(time/24, origin="1800-01-01")
  
  idx_lat <- which(lat >= box$lat[1] & lat <= box$lat[2])
  idx_lon <- c(which(lon >= box$lon[1]), which(lon <= box$lon[2]))
  
  # Carga DIRECTA sin buscar niveles (Start: Lon, Lat, Time)
  raw <- tryCatch({
    ncvar_get(nc, var, start=c(min(idx_lon), min(idx_lat), 1), 
              count=c(length(idx_lon), length(idx_lat), -1))
  }, error = function(e) return(NULL))
  
  if(is.null(raw)) return(NULL)
  
  # Promedio espacial
  series <- apply(raw, 3, mean, na.rm=TRUE)
  
  df <- data.frame(year=as.numeric(format(dates,"%Y")), month=as.numeric(format(dates,"%m")), val=series)
  df$clim_year <- ifelse(df$month==12, df$year+1, df$year)
  
  df %>% filter(month %in% c(12,1,2)) %>% 
    group_by(clim_year) %>% summarise(mean=mean(val, na.rm=TRUE)) %>% rename(year=clim_year)
}

# ==============================================================================
# 3. ROBUST DATA LOADING
# ==============================================================================
message("\n--- 1. CARGANDO DATOS (modo robusto) ---")
data_list <- list()

# 3.1 PRECIPITACIÓN (OBLIGATORIA)
if(file.exists(INPUT_PRECIP_FILE)) {
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
    dplyr::select(clim_year, PCP) %>%  
    rename(year=clim_year)
  
  data_list[["PCP"]] <- Y_dat
  message("✅ Precipitación: OK")
} else {
  stop("❌ ERROR: Falta archivo de precipitación.")
}

# 3.2 MOI (Mediterranean Oscillation)
if(file.exists(INPUT_MOI_FILE)) {
  raw_moi <- read.table(INPUT_MOI_FILE, header=TRUE)
  # Formato date: 18060101
  raw_moi$date_obj <- as.Date(as.character(raw_moi$date), format="%Y%m%d")
  raw_moi$year <- as.numeric(format(raw_moi$date_obj, "%Y"))
  raw_moi$month <- as.numeric(format(raw_moi$date_obj, "%m"))
  raw_moi$clim_year <- ifelse(raw_moi$month==12, raw_moi$year+1, raw_moi$year)
  
  X_moi <- raw_moi %>% filter(month %in% c(12,1,2)) %>%
    group_by(clim_year) %>% summarise(MOI = mean(moi, na.rm=TRUE)) %>% # 'moi' minúscula en tu archivo
    rename(year=clim_year)
  
  data_list[["MOI"]] <- X_moi
  message("✅ MOI: OK")
} else { message("⚠️ MOI: No encontrado") }

# 3.3 OTROS ÍNDICES (NAO, AMO, SOLAR, VOLCANES)
load_simple_index <- function(path, name, col_val) {
  if(file.exists(path)) {
    df <- read.table(path, header=TRUE)
    # Asumimos columnas estandarizadas o ajustamos aqui
    # Esta es una versión simplificada
    if("Year" %in% names(df)) df <- df %>% rename(year=Year)
    # Filtrar columnas
    df <- df %>% select(year, all_of(col_val))
    names(df)[2] <- name
    return(df)
  }
  return(NULL)
}

# Intentos de carga (Ajusta nombres de columnas según tus txt reales)
if(file.exists(INPUT_NAO_FILE)) data_list[["NAO"]] <- load_simple_index(INPUT_NAO_FILE, "NAO", "NAOmed")
if(file.exists(INPUT_SOLAR_FILE)) data_list[["SOLAR"]] <- load_simple_index(INPUT_SOLAR_FILE, "TSI", "TSI")
if(file.exists(INPUT_VOLC_FILE)) data_list[["VOLC"]] <- load_simple_index(INPUT_VOLC_FILE, "AOD", "AOD")

# 3.4 NETCDF (REANÁLISIS)
if(dir.exists(DIR_NETCDF)) {
  # U Wind Atlas
  tmp <- get_level(FILE_UWND, "uwnd", 300, BOX_ATLAS); if(!is.null(tmp)) data_list[["U_Atlas"]] <- tmp %>% rename(U_Atlas=mean)
  # U Wind Ireland (para Dipolo)
  tmp2 <- get_level(FILE_UWND, "uwnd", 300, BOX_IRELAND); 
  if(!is.null(tmp) && !is.null(tmp2)) {
    dip <- inner_join(tmp, tmp2, by="year") %>% mutate(DIP_Index = mean.x - mean.y) %>% select(year, DIP_Index)
    data_list[["DIP_Index"]] <- dip
  }
  # Pressure
  tmp <- get_sfc(FILE_PRMSL, "prmsl", BOX_ATLAS); if(!is.null(tmp)) data_list[["Press"]] <- tmp %>% rename(Press=mean)
  # RH 850
  tmp <- get_level(FILE_RHUM, "rhum", 850, BOX_ATLAS); if(!is.null(tmp)) data_list[["RH_850"]] <- tmp %>% rename(RH_850=mean)
  
  message("✅ NetCDF: Procesados los disponibles.")
} else { message("⚠️ NetCDF: Directorio no encontrado.") }

# 3.5 UNIÓN FINAL
dat <- Reduce(function(x, y) inner_join(x, y, by="year"), data_list)
dat <- na.omit(dat)
message(paste("\n>>> Dataset Final:", nrow(dat), "años. Variables:", paste(names(dat)[-1], collapse=", ")))

# ==============================================================================
# 4. DEFINICIÓN DE MODELOS (LINEAR, INTERACCIÓN, GAM)
# ==============================================================================

# Lista de fórmulas CANDIDATAS (Se filtrarán según datos disponibles)
candidate_formulas <- list(
  # --- UNIVARIATE ---
  "M01_Wind"        = "PCP ~ U_Atlas",
  "M02_Pres"        = "PCP ~ Press",
  "M03_Hum850"      = "PCP ~ RH_850",
  "M07_NAO"         = "PCP ~ NAO",
  "M09_MOI"         = "PCP ~ MOI",
  "M_Solar"         = "PCP ~ TSI",
  
  # --- MULTIVARIATE (ADDITIVE) ---
  "M10_Wind_Pres"   = "PCP ~ U_Atlas + Press",
  "M13_Sup_Base"    = "PCP ~ NAO + MOI",
  "M14_PRO_Local"   = "PCP ~ U_Atlas + Press + RH_850",
  "M19_ALL_STARS"   = "PCP ~ NAO + MOI + U_Atlas + Press + RH_850",
  
  # --- INTERACTIONS (SUPERVISOR TASK) ---
  # El ':' o '*' indica interacción. ¿Cambia el efecto de NAO según MOI?
  "M_Int_NAO_MOI"   = "PCP ~ NAO * MOI", 
  "M_Int_Wind_Pres" = "PCP ~ U_Atlas * Press",
  
  # --- GAMs (SUPERVISOR TASK: NO LINEAL / NO ESTACIONARIO) ---
  # s() indica "smooth term" (spline)
  "GAM_Trend"       = "PCP ~ s(year)",            # Tendencia temporal pura
  "GAM_Phy"         = "PCP ~ s(Press) + s(RH_850)", # Relaciones no lineales
  "GAM_Dyn"         = "PCP ~ s(NAO) + s(MOI)"
)

# Filtrado de fórmulas válidas
valid_models <- list()
avail_vars <- names(dat)

for(nm in names(candidate_formulas)) {
  f <- candidate_formulas[[nm]]
  # Extraer variables limpias
  vars <- all.vars(as.formula(f))
  vars <- vars[vars != "PCP" & vars != "s"] # Quitar Y y función spline
  
  if(all(vars %in% avail_vars)) {
    valid_models[[nm]] <- f
  }
}

# ==============================================================================
# 5. EJECUCIÓN Y STEPWISE
# ==============================================================================
results <- data.frame(Model=character(), Type=character(), R2_Adj=numeric(), AIC=numeric(), Equation=character(), stringsAsFactors=FALSE)

message("\n--- 2. EJECUTANDO MODELOS ---")

# A) MODELOS DEFINIDOS (LM y GAM)
for(name in names(valid_models)) {
  f_str <- valid_models[[name]]
  
  if(grepl("GAM", name)) {
    # --- Generalized Additive Model ---
    mod <- gam(as.formula(f_str), data=dat)
    s_mod <- summary(mod)
    
    r2  <- s_mod$r.sq
    aic <- AIC(mod)
    type <- "GAM"
    
    # --- CONSTRUIR "ECUACIÓN" GAM (Intercepto + EDF de cada variable) ---
    # 1. Intercepto
    beta0 <- coef(mod)[1]
    eq <- sprintf("%.2f", beta0)
    
    # 2. Términos Suavizados (Splines)
    if(!is.null(s_mod$s.table)) {
      # s_mod$s.table contiene los EDF en la primera columna
      for(k in 1:nrow(s_mod$s.table)) {
        var_name <- rownames(s_mod$s.table)[k] # Ej: "s(Press)"
        edf_val  <- s_mod$s.table[k, "edf"]    # El valor de curvatura
        
        # Formato: + s(Var)[edf=2.4]
        eq <- paste0(eq, sprintf(" + %s[edf=%.2f]", var_name, edf_val))
      }
    }
    
    # 3. Términos Lineales dentro del GAM (si los hubiera)
    if(!is.null(s_mod$p.table)) {
      # Saltamos el intercepto (row 1) si ya lo pusimos
      p_rows <- rownames(s_mod$p.table)
      p_rows <- p_rows[p_rows != "(Intercept)"]
      
      for(p_var in p_rows) {
        val <- s_mod$p.table[p_var, "Estimate"]
        signo <- ifelse(val >= 0, "+", "")
        eq <- paste0(eq, sprintf(" %s%.3f*%s", signo, val, p_var))
      }
    }
    
  } else {
    # --- Linear Model (Standard & Interaction) ---
    mod <- lm(as.formula(f_str), data=dat)
    r2  <- summary(mod)$adj.r.squared
    aic <- AIC(mod)
    type <- "LM"
    
    # Ecuación lineal estándar
    cf <- coef(mod)
    eq <- sprintf("%.2f", cf[1])
    for(i in 2:length(cf)) {
      if(is.na(cf[i])) next # Saltar coeficientes NA si los hubiera
      eq <- paste0(eq, sprintf(" %+0.3f*%s", cf[i], names(cf)[i]))
    }
  }
  
  results[nrow(results)+1,] <- list(name, type, r2*100, aic, eq)
}

# B) STEPWISE 
if(ncol(dat) > 3) {
  message(">>> Ejecutando Stepwise Selection...")
  full_formula <- as.formula(paste("PCP ~", paste(names(dat)[!names(dat) %in% c("year","PCP")], collapse=" + ")))
  full_mod <- lm(full_formula, data=dat)
  
  # Stepwise direction="both" busca el óptimo AIC
  step_mod <- stepAIC(full_mod, direction="both", trace=0)
  
  # --- EXTRAER ECUACIÓN (NUEVO) ---
  cf <- coef(step_mod)
  eq_step <- sprintf("%.2f", cf[1])
  for(i in 2:length(cf)) {
    val <- cf[i]
    signo <- ifelse(val >= 0, "+", "") # El negativo ya viene en el número
    var_name <- names(cf)[i]
    eq_step <- paste0(eq_step, sprintf(" %s%.4f*%s", signo, val, var_name))
  }
  
  results[nrow(results)+1,] <- list(
    "AUTO_Stepwise", 
    "STEP", 
    summary(step_mod)$adj.r.squared*100, 
    AIC(step_mod), 
    eq_step  
  )
}

# ==============================================================================
# 6. VISUALIZACIÓN Y RESULTADOS
# ==============================================================================
results$R2_Adj <- round(results$R2_Adj, 2)
results <- results %>% arrange(desc(R2_Adj))

# Imprimir Tabla
print(results %>% select(Model, Type, R2_Adj, AIC))

# Gráfico Lollipop Mejorado
results$Model <- factor(results$Model, levels=results$Model[order(results$R2_Adj)])

p <- ggplot(results, aes(x=R2_Adj, y=Model, color=Type)) +
  geom_segment(aes(x=0, xend=R2_Adj, y=Model, yend=Model), color="grey80", linetype="dotted") +
  geom_point(aes(size=2000/AIC), alpha=0.8) + # Tamaño inversamente prop. al AIC
  geom_text(aes(label=paste0(R2_Adj, "%")), hjust=-0.4, size=3) +
  scale_color_manual(values=c("LM"="steelblue", "GAM"="purple", "STEP"="orange")) +
  labs(title="Comparativa de modelos climáticos",
       subtitle="LM: Lineal | GAM: no lineal | STEP: automático",
       x="Varianza explicada (R2 adj %)", y="") +
  theme_minimal() +
  theme(legend.position="bottom")

print(p)

message("\n🏆 ECUACIÓN GANADORA (STEPWISE):")
print(results %>% filter(Type == "STEP") %>% pull(Equation))

message("\n✅ PROCESO FINALIZADO.")