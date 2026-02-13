# ==============================================================================
# 03_reconstruct_precip.R
# Author: Marcos Marin-Martin
# Date: 2026-02-11
# Description:
#   Advanced Dendroclimatic Reconstruction Pipeline.
#   
#   Methodology:
#   1. Input Flexibility: Handles Daily (cumulative window) or Monthly (seasonal sum) data.
#   2. Seasonality Logic: Automatically detects if the season crosses the calendar year
#      (e.g., Oct-Sep) or stays within the same year (e.g., May-Jul).
#   3. Dual Modeling: Performs both Linear (OLS) and Log-Normal (Scaling) reconstruction
#      to address heteroscedasticity and avoid negative precipitation values.
#   4. Validation: Calculates Split-period RE (Reduction of Error) and CE stats.
#
#   Outputs:
#   - _DATA.txt: Full series (Year, Std, Res, Obs, Rec_Lin, Rec_Log, MAs, SSS).
#   - _STATS.txt: Detailed model summaries and validation metrics.
# ==============================================================================

# --- 0. LOAD LIBRARIES ---
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(zoo)
  library(dplR)    
  library(ggplot2)
})

# ==============================================================================
# 1. CONFIGURATION & PARAMETERS
# ==============================================================================

# --- DATA TYPE ---
# 'daily'   = Uses daily sliding window (defined by WIN_DAYS)
# 'monthly' = Uses sum of specific months (defined by MONTH_START/END)
DATA_FREQ <- 'monthly' 

# --- MONTHLY PARAMETERS (Only if DATA_FREQ = 'monthly') ---
# Logic: 
# If START > END (e.g., 10 to 9), it calculates Hydrological Year (cross-year).
# If START <= END (e.g., 5 to 7), it calculates Calendar Season (same year).
MONTH_START <- 10  # Start Month
MONTH_END   <- 9   # End Month

# --- DAILY PARAMETERS (Only if DATA_FREQ = 'daily') ---
WIN_DAYS      <- 320
WIN_END_MONTH <- 6   
WIN_END_DAY   <- 30  

# --- FILE PATHS ---
DIR_BASE <- '/PLACEHOLDER/path/to/dir'

FILE_CLIMA <- file.path(DIR_BASE, 'PLACEHOLDER/path/to/climate.txt') 
FILE_CRONO <- file.path(DIR_BASE, 'PLACEHOLDER/path/to/chrono.txt')
FILE_RWL   <- file.path(DIR_BASE, 'PLACEHOLDER/path/to/chrono.rwi') 

# Output Files
OUT_DATA  <- file.path(DIR_BASE, 'Reconstructions/reconstruction.txt')
OUT_STATS <- file.path(DIR_BASE, 'Reconstructions/stats_models.txt')

# ==============================================================================
# 2. VALIDATION FUNCTIONS (RE/CE)
# ==============================================================================

calc_metrics <- function(obs, pred, cal_mean_obs, val_mean_obs) {
  # RE: Reduction of Error (vs calibration mean) - Critical for reconstruction validity
  RE <- 1 - sum((obs - pred)^2) / sum((obs - cal_mean_obs)^2)
  # CE: Coefficient of Efficiency (vs validation mean)
  CE <- 1 - sum((obs - pred)^2) / sum((obs - val_mean_obs)^2)
  return(c(RE=RE, CE=CE))
}

run_split_stats <- function(df, model_type="linear") {
  yrs <- sort(unique(df$year))
  mid <- median(yrs)
  early <- df %>% filter(year <= mid)
  late  <- df %>% filter(year > mid)
  
  # A) Calibrate Early -> Predict Late
  if(model_type=="linear") m_a <- lm(precip ~ ring, data=early) else m_a <- lm(log(precip) ~ ring, data=early)
  pred_late <- predict(m_a, newdata=late)
  if(model_type!="linear") pred_late <- exp(pred_late) 
  stats_late <- calc_metrics(late$precip, pred_late, mean(early$precip), mean(late$precip))
  
  # B) Calibrate Late -> Predict Early
  if(model_type=="linear") m_b <- lm(precip ~ ring, data=late) else m_b <- lm(log(precip) ~ ring, data=late)
  pred_early <- predict(m_b, newdata=early)
  if(model_type!="linear") pred_early <- exp(pred_early)
  stats_early <- calc_metrics(early$precip, pred_early, mean(late$precip), mean(early$precip))
  
  # Return average RE/CE
  return(list(
    RE = mean(c(stats_late["RE"], stats_early["RE"])),
    CE = mean(c(stats_late["CE"], stats_early["CE"]))
  ))
}

# ==============================================================================
# 3. CLIMATE DATA PROCESSING
# ==============================================================================

cat("\n>>> Processing Climate Data...\n")
precip_raw <- read.table(FILE_CLIMA, header = TRUE, sep = "\t")

# Ensure Date format
if("date" %in% names(precip_raw)) {
  precip_raw$Date <- as.Date(as.character(precip_raw$date), format="%Y%m%d")
  precip_raw$Year <- year(precip_raw$Date)
  precip_raw$Month <- month(precip_raw$Date)
} else if ("Year" %in% names(precip_raw) & "Month" %in% names(precip_raw)) {
  # Format is already monthly columns
} else {
  stop("Input file must have 'date' column or 'Year'/'Month' columns.")
}

clim_annual <- data.frame()

if (DATA_FREQ == 'daily') {
  # Daily logic (Sliding Window)
  cat(paste(">>> Mode: DAILY | Window:", WIN_DAYS, "days.\n"))
  min_dt <- min(precip_raw$Date); max_dt <- max(precip_raw$Date)
  start_yr <- year(min_dt) + 1; end_yr <- year(max_dt)
  
  res_list <- list()
  for(y in start_yr:end_yr) {
    end_date <- ymd(sprintf("%d-%02d-%02d", y, WIN_END_MONTH, WIN_END_DAY))
    start_date <- end_date - days(WIN_DAYS - 1)
    if(start_date >= min_dt && end_date <= max_dt) {
      subset_p <- precip_raw %>% filter(Date >= start_date & Date <= end_date)
      if(nrow(subset_p) >= (WIN_DAYS - 5)) {
        res_list[[as.character(y)]] <- data.frame(year=y, precip=sum(subset_p$precipitation, na.rm=TRUE))
      }
    }
  }
  clim_annual <- bind_rows(res_list)
  
} else if (DATA_FREQ == 'monthly') {
  # Monthly logic (Seasonal Sum)
  is_cross_year <- MONTH_START > MONTH_END
  cat(paste(">>> Mode: MONTHLY | Season:", MONTH_START, "to", MONTH_END, 
            "| Cross-Year:", is_cross_year, "\n"))
  
  clim_annual <- precip_raw %>%
    mutate(
      # Assign Hydrological Year if crossing Dec (e.g., Oct 1950 -> HydroYear 1951)
      HydroYear = if(is_cross_year) {
        ifelse(Month >= MONTH_START, Year + 1, Year)
      } else {
        Year
      }
    ) %>%
    filter(
      if (is_cross_year) {
        (Month >= MONTH_START) | (Month <= MONTH_END)
      } else {
        (Month >= MONTH_START) & (Month <= MONTH_END)
      }
    ) %>%
    group_by(HydroYear) %>%
    summarise(precip = sum(precipitation, na.rm=TRUE), n_months = n()) %>%
    # Determine expected number of months
    filter(n_months == (ifelse(is_cross_year, 12 - MONTH_START + 1 + MONTH_END, MONTH_END - MONTH_START + 1))) %>%
    select(year = HydroYear, precip)
}

cat(paste(">>> Climate Years Available:", nrow(clim_annual), "\n"))

# ==============================================================================
# 4. LOAD CHRONOLOGIES & SSS
# ==============================================================================

crono_df <- read.table(FILE_CRONO, header = TRUE, sep = "\t")

cat(">>> Calculating SSS from RWL...\n")
sss_vec <- rep(NA, nrow(crono_df))
if(file.exists(FILE_RWL)) {
  tryCatch({
    rwl <- read.rwl(FILE_RWL)
    sss_val <- sss(rwl)
    idx <- match(crono_df$year, as.numeric(names(sss_val)))
    sss_vec <- sss_val[idx]
  }, error=function(e) cat("Warning: Could not calculate SSS.\n"))
}
crono_df$SSS <- sss_vec

# Join for Calibration
df_calib <- inner_join(clim_annual, crono_df, by="year")

if(nrow(df_calib) < 10) stop("CRITICAL ERROR: Less than 10 overlapping years for calibration.")

# ==============================================================================
# 5. MODELING & STATISTICS REPORT
# ==============================================================================

stats_text <- c(
  "============================================================",
  paste(" STATISTICAL REPORT: DENDROCLIMATIC RECONSTRUCTION"),
  paste(" Date:", Sys.Date()),
  paste(" Frequency:", DATA_FREQ),
  "============================================================"
)

# We will save the models based on 'res' chronology for the final output
final_mod_lin <- NULL
final_mod_log <- NULL

for (tipo in c("std", "res")) {
  
  curr_data <- df_calib %>% select(year, precip, ring = all_of(tipo))
  
  stats_text <- c(stats_text, paste("\n\n################ CHRONOLOGY TYPE:", toupper(tipo), "################"))
  
  # --- A) LINEAR MODEL (OLS) ---
  mod_lin <- lm(precip ~ ring, data = curr_data)
  s_lin   <- summary(mod_lin)
  val_lin <- run_split_stats(curr_data, "linear")
  
  stats_text <- c(stats_text, 
                  "\n--- [1] LINEAR MODEL (OLS) ---",
                  paste("Equation: P =", round(coef(mod_lin)[1],4), "+", round(coef(mod_lin)[2],4), "* Index"),
                  paste("RMSE (Cal):", round(s_lin$sigma, 2)),
                  paste("R2:", round(s_lin$r.squared, 4), "| Adj-R2:", round(s_lin$adj.r.squared, 4)),
                  paste("RE (Split):", round(val_lin$RE, 4), "| CE (Split):", round(val_lin$CE, 4)),
                  "\n> R Summary (Linear):")
  stats_text <- c(stats_text, capture.output(s_lin))
  
  # --- B) LOG-NORMAL MODEL (Scaling) ---
  # log(P) = a + b * Index
  mod_log <- lm(log(precip) ~ ring, data = curr_data)
  s_log   <- summary(mod_log)
  
  # Calculate Stats in Real Space (mm) for valid comparison
  pred_mm <- exp(predict(mod_log))
  obs_mm  <- curr_data$precip
  r2_real <- cor(obs_mm, pred_mm)^2
  rmse_real <- sqrt(mean((obs_mm - pred_mm)^2))
  val_log <- run_split_stats(curr_data, "log")
  
  stats_text <- c(stats_text, 
                  "\n--- [2] LOG-NORMAL MODEL (Log-Scaling) ---",
                  paste("Equation: P = exp(", round(coef(mod_log)[1],4), "+", round(coef(mod_log)[2],4), "* Index )"),
                  "NOTE: R2/RMSE below are back-transformed to real space (mm).",
                  paste("R2 (Real mm):", round(r2_real, 4)),
                  paste("RMSE (Real mm):", round(rmse_real, 2)),
                  paste("RE (Split):", round(val_log$RE, 4), "| CE (Split):", round(val_log$CE, 4)),
                  "\n> R Summary (Log-Space):")
  stats_text <- c(stats_text, capture.output(s_log))
  
  # Save 'res' models for final reconstruction
  if(tipo == "res") {
    final_mod_lin <- mod_lin
    final_mod_log <- mod_log
  }
}

# Write Stats to File
writeLines(stats_text, OUT_STATS)
cat(paste(">>> Statistics saved to:", OUT_STATS, "\n"))

# ==============================================================================
# 6. FINAL RECONSTRUCTION OUTPUT
# ==============================================================================

# Using 'res' chronology for the final output file
rec_df <- crono_df %>%
  mutate(
    # Merge observed precipitation
    precipitation = clim_annual$precip[match(year, clim_annual$year)], 
    
    # 1. Linear Reconstruction
    rec_lineal = predict(final_mod_lin, newdata = data.frame(ring = res)),
    
    # 2. Log Reconstruction (Back-transform)
    rec_log = exp(predict(final_mod_log, newdata = data.frame(ring = res))),
    
    # 3. Moving Averages (11-yr)
    `11y_ma_lineal` = rollmean(rec_lineal, k=11, fill=NA, align="center"),
    `11y_ma_log`    = rollmean(rec_log, k=11, fill=NA, align="center"),
    
    # 4. Moving Averages (31-yr)
    `31y_ma_lineal` = rollmean(rec_lineal, k=31, fill=NA, align="center"),
    `31y_ma_log`    = rollmean(rec_log, k=31, fill=NA, align="center")
  ) %>%
  # Select and Order Columns
  select(year, std, res, precipitation, 
         rec_lineal, rec_log, 
         `11y_ma_lineal`, `11y_ma_log`, 
         `31y_ma_lineal`, `31y_ma_log`, 
         SSS)

# Write Data to File
write.table(rec_df, OUT_DATA, sep="\t", row.names=FALSE, quote=FALSE, na="NA")
cat(paste(">>> Full Reconstruction Data saved to:", OUT_DATA, "\n"))

# ==============================================================================
# 7. QUICK PLOT
# ==============================================================================
p <- ggplot(rec_df, aes(x=year)) +
  geom_line(aes(y=rec_log, color="Rec (Log-Normal)"), size=0.4, alpha=0.7) +
  geom_line(aes(y=`11y_ma_log`, color="MA 11-yr"), size=0.9) +
  geom_line(aes(y=precipitation, color="Observed"), size=0.6) +
  scale_color_manual(values=c("Rec (Log-Normal)"="steelblue", "MA 11-yr"="red", "Observed"="black")) +
  labs(title="Final Reconstruction", subtitle=paste("Frequency:", DATA_FREQ), y="Precipitation (mm)") +
  theme_minimal() + theme(legend.position="bottom")

print(p)
cat(">>> SCRIPT COMPLETED SUCCESSFULLY.\n")