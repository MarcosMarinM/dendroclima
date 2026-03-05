# ==============================================================================
# 01_main_reconstruction.R
# Author: Marcos Marín-Martín
# Date: 2026-03-05
# Description:
#   Advanced Dendroclimatic Reconstruction Pipeline.
#
#   Methodology:
#   1. Input Flexibility: Handles Daily (cumulative window) or Monthly (seasonal
#      sum) data. Automatically detects cross-year seasons (e.g., Oct-Sep).
#   2. Dual Modeling: Performs both Linear (OLS) and Log-Normal reconstruction
#      to address heteroscedasticity and avoid negative precipitation values.
#   3. Validation: Calculates split-period RE and CE statistics.
#   4. Quality: SSS (Subsample Signal Strength) calculated and used to filter
#      the reconstruction plot below a configurable threshold.
#   5. Visualisation: 11-yr smoothing, RMSE uncertainty ribbon, horizontal mean
#      line, and per-century extreme year labels.
#
#   Inputs:
#   - Monthly or daily climate data file (.txt, tab-separated).
#   - Residual chronology file (.txt, tab-separated).
#   - RWI file for SSS calculation (.rwi).
#
#   Outputs:
#   - _DATA.txt     : Full reconstruction series (year, std, res, precipitation,
#                     rec_lineal, rec_log, moving averages, SSS).
#   - _STATS.txt    : Detailed model summaries and split-period validation.
#   - _EQUATION.txt : Model equations, R², RMSE.
#   - _PLOT.pdf     : Reconstruction plot with ribbon, smoothing, extremes.
# ==============================================================================

# --- 0. LIBRARIES ---
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(zoo)
  library(dplR)
  library(ggplot2)
  library(ggrepel)   # for non-overlapping extreme year labels (install if needed)
})

# ==============================================================================
# 1. CONFIGURATION & PARAMETERS
# ==============================================================================

# --- DATA TYPE ---
# 'daily'   = Uses daily sliding window (defined by WIN_DAYS)
# 'monthly' = Uses sum of specific months (defined by MONTH_START/END)
DATA_FREQ <- 'monthly'

# --- MONTHLY PARAMETERS (only if DATA_FREQ = 'monthly') ---
# If START > END (e.g., 10 to 9)  -> Hydrological Year (cross-year).
# If START <= END (e.g., 5 to 7)  -> Calendar Season (same year).
MONTH_START <- 10
MONTH_END   <- 9

# --- DAILY PARAMETERS (only if DATA_FREQ = 'daily') ---
WIN_DAYS      <- 320
WIN_END_MONTH <- 6
WIN_END_DAY   <- 30

# --- PLOT OPTIONS ---
# Minimum SSS to include a year in the reconstruction plot (set to 0 to disable)
SSS_THRESHOLD <- 0.85
# Which model to use for the main plot: "linear" or "log"
PLOT_MODEL    <- "log"

# --- FILE PATHS ---
DIR_BASE <- 'PLACEHOLDER/path/to/dir'

FILE_CLIMA <- file.path(DIR_BASE, 'PLACEHOLDER/path/to/climate.txt')
FILE_CRONO <- file.path(DIR_BASE, 'PLACEHOLDER/path/to/chrono.txt')
FILE_RWL   <- file.path(DIR_BASE, 'PLACEHOLDER/path/to/chrono.rwi')

OUT_DATA     <- file.path(DIR_BASE, 'Reconstructions/reconstruction_DATA.txt')
OUT_STATS    <- file.path(DIR_BASE, 'Reconstructions/reconstruction_STATS.txt')
OUT_EQUATION <- file.path(DIR_BASE, 'Reconstructions/reconstruction_EQUATION.txt')
OUT_PLOT     <- file.path(DIR_BASE, 'Reconstructions/reconstruction_PLOT.pdf')

if (!dir.exists(dirname(OUT_DATA))) dir.create(dirname(OUT_DATA), recursive = TRUE)

# ==============================================================================
# 2. VALIDATION FUNCTIONS (RE / CE)
# ==============================================================================

calc_metrics <- function(obs, pred, cal_mean_obs, val_mean_obs) {
  RE <- 1 - sum((obs - pred)^2) / sum((obs - cal_mean_obs)^2)
  CE <- 1 - sum((obs - pred)^2) / sum((obs - val_mean_obs)^2)
  c(RE = RE, CE = CE)
}

run_split_stats <- function(df, model_type = "linear") {
  yrs <- sort(unique(df$year))
  mid <- median(yrs)
  early <- df |> filter(year <= mid)
  late  <- df |> filter(year >  mid)

  fit_model <- function(train, test) {
    if (model_type == "linear") {
      m <- lm(precip ~ ring, data = train)
      pred <- predict(m, newdata = test)
    } else {
      m <- lm(log(precip) ~ ring, data = train)
      pred <- exp(predict(m, newdata = test))
    }
    list(m = m, pred = pred)
  }

  res_a <- fit_model(early, late)
  res_b <- fit_model(late, early)

  stats_a <- calc_metrics(late$precip,  res_a$pred, mean(early$precip), mean(late$precip))
  stats_b <- calc_metrics(early$precip, res_b$pred, mean(late$precip),  mean(early$precip))

  list(
    RE = mean(c(stats_a["RE"], stats_b["RE"])),
    CE = mean(c(stats_a["CE"], stats_b["CE"]))
  )
}

# ==============================================================================
# 3. CLIMATE DATA PROCESSING
# ==============================================================================

cat("\n>>> Processing Climate Data...\n")
precip_raw <- read.table(FILE_CLIMA, header = TRUE, sep = "\t")

if ("date" %in% names(precip_raw)) {
  precip_raw$Date  <- as.Date(as.character(precip_raw$date), format = "%Y%m%d")
  precip_raw$Year  <- year(precip_raw$Date)
  precip_raw$Month <- month(precip_raw$Date)
} else if ("Year" %in% names(precip_raw) & "Month" %in% names(precip_raw)) {
  # Already in monthly format
} else {
  stop("Input climate file must have a 'date' column (YYYYMMDD) or 'Year'/'Month' columns.")
}

clim_annual <- data.frame()

if (DATA_FREQ == 'daily') {
  cat(paste(">>> Mode: DAILY | Window:", WIN_DAYS, "days.\n"))
  min_dt <- min(precip_raw$Date); max_dt <- max(precip_raw$Date)
  start_yr <- year(min_dt) + 1;   end_yr   <- year(max_dt)

  res_list <- list()
  for (y in start_yr:end_yr) {
    end_date   <- ymd(sprintf("%d-%02d-%02d", y, WIN_END_MONTH, WIN_END_DAY))
    start_date <- end_date - days(WIN_DAYS - 1)
    if (start_date >= min_dt && end_date <= max_dt) {
      sub <- precip_raw |> filter(Date >= start_date & Date <= end_date)
      if (nrow(sub) >= (WIN_DAYS - 5)) {
        res_list[[as.character(y)]] <- data.frame(year = y,
                                                   precip = sum(sub$precipitation, na.rm = TRUE))
      }
    }
  }
  clim_annual <- bind_rows(res_list)

} else if (DATA_FREQ == 'monthly') {
  is_cross_year <- MONTH_START > MONTH_END
  cat(paste(">>> Mode: MONTHLY | Season:", MONTH_START, "to", MONTH_END,
            "| Cross-Year:", is_cross_year, "\n"))

  expected_months <- if (is_cross_year) {
    12 - MONTH_START + 1 + MONTH_END
  } else {
    MONTH_END - MONTH_START + 1
  }

  clim_annual <- precip_raw |>
    mutate(
      HydroYear = if (is_cross_year) ifelse(Month >= MONTH_START, Year + 1, Year) else Year
    ) |>
    filter(
      if (is_cross_year) (Month >= MONTH_START) | (Month <= MONTH_END)
      else               (Month >= MONTH_START) & (Month <= MONTH_END)
    ) |>
    group_by(HydroYear) |>
    summarise(precip = sum(precipitation, na.rm = TRUE), n_months = n(), .groups = "drop") |>
    filter(n_months == expected_months) |>
    select(year = HydroYear, precip)
}

cat(paste(">>> Climate years available:", nrow(clim_annual), "\n"))

# ==============================================================================
# 4. LOAD CHRONOLOGIES & SSS
# ==============================================================================

crono_df <- read.table(FILE_CRONO, header = TRUE, sep = "\t")

# Normalise year column name
if (!"year" %in% names(crono_df) && "Year" %in% names(crono_df)) {
  crono_df <- crono_df |> rename(year = Year)
}
if (!"res" %in% names(crono_df)) stop("Chronology file must have a 'res' column.")
if (!"std" %in% names(crono_df)) stop("Chronology file must have a 'std' column.")

crono_df$year <- as.numeric(crono_df$year)

cat(">>> Calculating SSS from RWL...\n")
sss_vec <- rep(NA_real_, nrow(crono_df))
if (file.exists(FILE_RWL)) {
  tryCatch({
    rwl     <- read.rwl(FILE_RWL)
    sss_val <- sss(rwl)
    idx     <- match(crono_df$year, as.numeric(names(sss_val)))
    sss_vec <- sss_val[idx]
  }, error = function(e) cat("Warning: Could not calculate SSS:", e$message, "\n"))
}
crono_df$SSS <- sss_vec

# Join for calibration
df_calib <- inner_join(clim_annual, crono_df, by = "year")
if (nrow(df_calib) < 10) stop("Less than 10 overlapping years. Cannot calibrate.")
cat(paste(">>> Calibration overlap:", nrow(df_calib), "years.\n"))

# ==============================================================================
# 5. MODELLING & STATISTICS REPORT
# ==============================================================================

stats_text <- c(
  "============================================================",
  "  STATISTICAL REPORT: DENDROCLIMATIC RECONSTRUCTION",
  paste("  Date:", Sys.Date()),
  paste("  Frequency:", DATA_FREQ),
  "============================================================"
)

final_mod_lin <- NULL
final_mod_log <- NULL

for (tipo in c("std", "res")) {
  curr_data <- df_calib |> select(year, precip, ring = all_of(tipo))

  stats_text <- c(stats_text,
                  paste("\n\n################ CHRONOLOGY TYPE:", toupper(tipo), "################"))

  # A) Linear OLS
  mod_lin <- lm(precip ~ ring, data = curr_data)
  s_lin   <- summary(mod_lin)
  val_lin <- run_split_stats(curr_data, "linear")

  stats_text <- c(stats_text,
    "\n--- [1] LINEAR MODEL (OLS) ---",
    paste("Equation: P =", round(coef(mod_lin)[1], 4), "+", round(coef(mod_lin)[2], 4), "* Index"),
    paste("RMSE (Cal):", round(s_lin$sigma, 2)),
    paste("R2:", round(s_lin$r.squared, 4), "| Adj-R2:", round(s_lin$adj.r.squared, 4)),
    paste("RE (Split):", round(val_lin$RE, 4), "| CE (Split):", round(val_lin$CE, 4)),
    "\n> R Summary (Linear):",
    capture.output(s_lin)
  )

  # B) Log-Normal
  mod_log   <- lm(log(precip) ~ ring, data = curr_data)
  s_log     <- summary(mod_log)
  pred_mm   <- exp(predict(mod_log))
  obs_mm    <- curr_data$precip
  r2_real   <- cor(obs_mm, pred_mm)^2
  rmse_real <- sqrt(mean((obs_mm - pred_mm)^2))
  val_log   <- run_split_stats(curr_data, "log")

  stats_text <- c(stats_text,
    "\n--- [2] LOG-NORMAL MODEL ---",
    paste("Equation: P = exp(", round(coef(mod_log)[1], 4), "+",
          round(coef(mod_log)[2], 4), "* Index )"),
    "NOTE: R2/RMSE back-transformed to real space (mm).",
    paste("R2 (Real mm):", round(r2_real, 4)),
    paste("RMSE (Real mm):", round(rmse_real, 2)),
    paste("RE (Split):", round(val_log$RE, 4), "| CE (Split):", round(val_log$CE, 4)),
    "\n> R Summary (Log-Space):",
    capture.output(s_log)
  )

  if (tipo == "res") {
    final_mod_lin <- mod_lin
    final_mod_log <- mod_log
  }
}

writeLines(stats_text, OUT_STATS)
cat(paste(">>> Statistics saved to:", OUT_STATS, "\n"))

# ==============================================================================
# 6. EQUATION FILE
# ==============================================================================

rmse_lin <- sqrt(mean((df_calib$precip - predict(final_mod_lin))^2))
rmse_log <- sqrt(mean((df_calib$precip - exp(predict(final_mod_log)))^2))

eq_lines <- c(
  paste("Date:", Sys.Date()),
  paste("Chronology type: res | Frequency:", DATA_FREQ),
  "",
  "--- LINEAR MODEL ---",
  sprintf("P = %.4f + %.4f * res", coef(final_mod_lin)[1], coef(final_mod_lin)[2]),
  sprintf("R2   = %.4f", summary(final_mod_lin)$r.squared),
  sprintf("RMSE = %.2f mm", rmse_lin),
  "",
  "--- LOG-NORMAL MODEL ---",
  sprintf("P = exp( %.4f + %.4f * res )", coef(final_mod_log)[1], coef(final_mod_log)[2]),
  sprintf("R2   = %.4f  (back-transformed to mm)", cor(df_calib$precip, exp(predict(final_mod_log)))^2),
  sprintf("RMSE = %.2f mm  (back-transformed)", rmse_log)
)

writeLines(eq_lines, OUT_EQUATION)
cat(paste(">>> Equation file saved to:", OUT_EQUATION, "\n"))

# ==============================================================================
# 7. FULL RECONSTRUCTION OUTPUT
# ==============================================================================

rec_df <- crono_df |>
  mutate(
    precipitation   = clim_annual$precip[match(year, clim_annual$year)],
    rec_lineal      = predict(final_mod_lin, newdata = data.frame(ring = res)),
    rec_log         = exp(predict(final_mod_log, newdata = data.frame(ring = res))),
    `11y_ma_lineal` = rollmean(rec_lineal, k = 11, fill = NA, align = "center"),
    `11y_ma_log`    = rollmean(rec_log,    k = 11, fill = NA, align = "center"),
    `31y_ma_lineal` = rollmean(rec_lineal, k = 31, fill = NA, align = "center"),
    `31y_ma_log`    = rollmean(rec_log,    k = 31, fill = NA, align = "center")
  ) |>
  select(year, std, res, precipitation,
         rec_lineal, rec_log,
         `11y_ma_lineal`, `11y_ma_log`,
         `31y_ma_lineal`, `31y_ma_log`,
         SSS)

write.table(rec_df, OUT_DATA, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
cat(paste(">>> Full reconstruction saved to:", OUT_DATA, "\n"))

# ==============================================================================
# 8. PLOT
# ==============================================================================

rec_col    <- if (PLOT_MODEL == "log") "rec_log"         else "rec_lineal"
ma_col     <- if (PLOT_MODEL == "log") "11y_ma_log"      else "11y_ma_lineal"
rmse_plot  <- if (PLOT_MODEL == "log") rmse_log          else rmse_lin

# Filter by SSS threshold
plot_df <- rec_df |>
  filter(is.na(SSS) | SSS >= SSS_THRESHOLD) |>
  rename(rec_plot = all_of(rec_col), ma_plot = all_of(ma_col))

# Per-century extremes (labels)
plot_df <- plot_df |> mutate(Century = (year %/% 100) * 100)

low_yrs  <- plot_df |> group_by(Century) |>
  filter(rec_plot == min(rec_plot, na.rm = TRUE)) |> slice(1)
high_yrs <- plot_df |> group_by(Century) |>
  filter(rec_plot == max(rec_plot, na.rm = TRUE)) |> slice(1)

mean_obs <- mean(plot_df$precipitation, na.rm = TRUE)

p <- ggplot(plot_df, aes(x = year)) +
  geom_ribbon(aes(ymin = rec_plot - rmse_plot, ymax = rec_plot + rmse_plot),
              fill = "grey70", alpha = 0.4) +
  geom_line(aes(y = rec_plot,     colour = "Reconstructed"), linewidth = 0.4, alpha = 0.8) +
  geom_line(aes(y = precipitation, colour = "Observed"),    linewidth = 0.5, linetype = "dashed") +
  geom_line(aes(y = ma_plot),      colour = "#8B6969",      linewidth = 0.7) +
  geom_hline(yintercept = mean_obs, linetype = "dashed", colour = "#CDAD00", linewidth = 0.5) +
  geom_text(data = low_yrs,  aes(x = year, y = rec_plot, label = year),
            angle = 45, vjust = 1.8, hjust = 1,  size = 2.5, colour = "red") +
  geom_text(data = high_yrs, aes(x = year, y = rec_plot, label = year),
            angle = 45, vjust = -0.8, hjust = 0, size = 2.5, colour = "blue") +
  scale_colour_manual(values = c("Reconstructed" = "steelblue", "Observed" = "black")) +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  labs(
    x      = if (DATA_FREQ == "monthly" && MONTH_START > MONTH_END)
                sprintf("Hydrological year (%s–%s)", month.abb[MONTH_START], month.abb[MONTH_END])
              else "Year",
    y      = "Precipitation (mm)",
    colour = NULL,
    caption = sprintf("Grey ribbon = ±RMSE (%.1f mm) | Brown line = 11-yr MA | SSS ≥ %.2f",
                      rmse_plot, SSS_THRESHOLD)
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(p, filename = OUT_PLOT, width = 12, height = 5)
cat(paste(">>> Plot saved to:", OUT_PLOT, "\n"))

cat("\n>>> SCRIPT COMPLETED SUCCESSFULLY.\n")
cat("    Files generated:\n")
cat("    -", OUT_DATA,     "\n")
cat("    -", OUT_STATS,    "\n")
cat("    -", OUT_EQUATION, "\n")
cat("    -", OUT_PLOT,     "\n")
