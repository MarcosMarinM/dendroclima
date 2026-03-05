# ==============================================================================
# 02_calibration_verification.R
# Author: Marcos Marín-Martín
# Date: 2026-03-05
# Description:
#   Split-Period Calibration & Verification Analysis.
#   Validates the stability and predictive skill of the reconstruction model by
#   splitting the instrumental period into two equal halves (Early / Late).
#
#   Procedure:
#   - Round A: Calibrate on Early Period -> Verify on Late Period.
#   - Round B: Calibrate on Late Period -> Verify on Early Period.
#
#   Metrics:
#   - R² (calibration and verification)
#   - Reduction of Error (RE)
#   - Coefficient of Efficiency (CE)
#
#   Input:
#   - _DATA.txt file produced by 03_reconstruct_precip.R
#     Required columns: year, res, std, precipitation, rec_lineal, rec_log
#
#   Output:
#   - 4-panel PDF plot (Calibration/Verification for both rounds)
# ==============================================================================

# --- 0. LIBRARIES ---
suppressPackageStartupMessages({
  library(dplyr)
})

# ==============================================================================
# 1. PARAMETERS
# ==============================================================================

INPUT_FILE_RECON <- 'PLACEHOLDER/path/to/reconstruction_DATA.txt'
OUTPUT_PDF_FILE  <- 'PLACEHOLDER/path/to/calibration_verification_plot.pdf'

# Which chronology index to use as predictor: "res" or "std"
CHRONO_TYPE <- "res"

# Which model to validate: "linear" or "log"
MODEL_TYPE <- "linear"

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

calc_r2 <- function(obs, pred) {
  cor(obs, pred, use = "complete.obs")^2
}

calc_re <- function(obs, pred, cal_mean) {
  # Reduction of Error: skill vs. calibration-period mean as benchmark
  1 - sum((obs - pred)^2) / sum((obs - cal_mean)^2)
}

calc_ce <- function(obs, pred, ver_mean) {
  # Coefficient of Efficiency: skill vs. verification-period mean as benchmark
  1 - sum((obs - pred)^2) / sum((obs - ver_mean)^2)
}

fit_and_predict <- function(cal_df, ver_df, model_type) {
  if (model_type == "linear") {
    m <- lm(precipitation ~ ring, data = cal_df)
    pred <- predict(m, newdata = ver_df)
  } else {
    m <- lm(log(precipitation) ~ ring, data = cal_df)
    pred <- exp(predict(m, newdata = ver_df))
  }
  list(model = m, pred = pred)
}

# ==============================================================================
# 3. LOAD DATA
# ==============================================================================

tryCatch({
  raw <- read.table(INPUT_FILE_RECON, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)
  cat("Data loaded correctly. Columns found:", paste(names(raw), collapse = ", "), "\n")
}, error = function(e) {
  stop("Error loading file. Check path and separator. Error: ", e$message)
})

# Coerce types
raw$year          <- as.numeric(raw$year)
raw$precipitation <- as.numeric(raw$precipitation)
raw$res           <- as.numeric(raw$res)
raw$std           <- as.numeric(raw$std)

# ==============================================================================
# 4. PREPARE CALIBRATION DATASET
# ==============================================================================

# Keep only the instrumental period (years with observed precipitation)
df_inst <- raw |>
  filter(!is.na(precipitation)) |>
  select(year, precipitation, ring = all_of(CHRONO_TYPE)) |>
  filter(!is.na(ring)) |>
  arrange(year)

if (nrow(df_inst) < 20) {
  stop("Less than 20 overlapping years found. Cannot perform split-period validation.")
}

# Split at median year (equal halves)
mid_year   <- median(df_inst$year)
period1    <- df_inst |> filter(year <= mid_year)   # Early
period2    <- df_inst |> filter(year >  mid_year)   # Late

period1_years <- range(period1$year)
period2_years <- range(period2$year)

cat(sprintf("Instrumental period: %d – %d (%d years)\n",
            min(df_inst$year), max(df_inst$year), nrow(df_inst)))
cat(sprintf("Early period: %d – %d (%d years)\n",
            period1_years[1], period1_years[2], nrow(period1)))
cat(sprintf("Late period:  %d – %d (%d years)\n",
            period2_years[1], period2_years[2], nrow(period2)))

# ==============================================================================
# 5. ROUND A: Calibrate Early -> Verify Late
# ==============================================================================

fit_a   <- fit_and_predict(period1, period2, MODEL_TYPE)
pred_a  <- fit_and_predict(period1, period1, MODEL_TYPE)  # calibration fit on itself

r2_cal1 <- calc_r2(period1$precipitation, pred_a$pred)
r2_ver1 <- calc_r2(period2$precipitation, fit_a$pred)
re_ver1 <- calc_re(period2$precipitation, fit_a$pred, mean(period1$precipitation))
ce_ver1 <- calc_ce(period2$precipitation, fit_a$pred, mean(period2$precipitation))

results_cal1 <- data.frame(
  year              = period1$year,
  observed          = period1$precipitation,
  predicted         = pred_a$pred
)

results_ver1 <- data.frame(
  year              = period2$year,
  observed          = period2$precipitation,
  predicted         = fit_a$pred
)

# ==============================================================================
# 6. ROUND B: Calibrate Late -> Verify Early
# ==============================================================================

fit_b   <- fit_and_predict(period2, period1, MODEL_TYPE)
pred_b  <- fit_and_predict(period2, period2, MODEL_TYPE)  # calibration fit on itself

r2_cal2 <- calc_r2(period2$precipitation, pred_b$pred)
r2_ver2 <- calc_r2(period1$precipitation, fit_b$pred)
re_ver2 <- calc_re(period1$precipitation, fit_b$pred, mean(period2$precipitation))
ce_ver2 <- calc_ce(period1$precipitation, fit_b$pred, mean(period1$precipitation))

results_cal2 <- data.frame(
  year              = period2$year,
  observed          = period2$precipitation,
  predicted         = pred_b$pred
)

results_ver2 <- data.frame(
  year              = period1$year,
  observed          = period1$precipitation,
  predicted         = fit_b$pred
)

# ==============================================================================
# 7. PRINT SUMMARY
# ==============================================================================

cat("\n========================================\n")
cat("  SPLIT-PERIOD CALIBRATION / VERIFICATION\n")
cat("========================================\n")
cat(sprintf("Model type   : %s\n", toupper(MODEL_TYPE)))
cat(sprintf("Chrono index : %s\n\n", toupper(CHRONO_TYPE)))

cat(sprintf("ROUND A  |  Calibration (%d–%d)  ->  Verification (%d–%d)\n",
            period1_years[1], period1_years[2], period2_years[1], period2_years[2]))
cat(sprintf("  R² (cal) = %.3f   R² (ver) = %.3f\n", r2_cal1, r2_ver1))
cat(sprintf("  RE (ver) = %.3f   CE (ver) = %.3f\n\n", re_ver1, ce_ver1))

cat(sprintf("ROUND B  |  Calibration (%d–%d)  ->  Verification (%d–%d)\n",
            period2_years[1], period2_years[2], period1_years[1], period1_years[2]))
cat(sprintf("  R² (cal) = %.3f   R² (ver) = %.3f\n", r2_cal2, r2_ver2))
cat(sprintf("  RE (ver) = %.3f   CE (ver) = %.3f\n", re_ver2, ce_ver2))
cat("========================================\n")

# ==============================================================================
# 8. PLOT
# ==============================================================================

# Shared y-axis range across all panels
all_vals <- c(results_cal1$observed, results_cal1$predicted,
              results_ver1$observed, results_ver1$predicted,
              results_cal2$observed, results_cal2$predicted,
              results_ver2$observed, results_ver2$predicted)
y_range <- range(all_vals, na.rm = TRUE)
y_range <- y_range + c(-1, 1) * diff(y_range) * 0.05  # 5% padding

base_pointsize <- 10
cex_axis_val   <- 0.8
cex_lab_val    <- 0.9
cex_title_val  <- 0.75
cex_legend_val <- 0.8

format_r2 <- function(x) if (is.na(x)) "NA" else sprintf("%.2f", x)
format_rece <- function(re, ce) sprintf("RE=%.2f, CE=%.2f", re, ce)

pdf(OUTPUT_PDF_FILE, width = 8, height = 5,
    pointsize = base_pointsize, useDingbats = FALSE)

par(mfrow = c(2, 2),
    mar   = c(3.5, 3.5, 1.5, 1),
    oma   = c(0, 0, 0, 0),
    mgp   = c(1.8, 0.6, 0),
    cex.axis = cex_axis_val,
    cex.lab  = cex_lab_val)

plot_panel <- function(df, title_main, subtitle, col_pred = "red") {
  plot(df$year, df$observed, type = "l", col = "black",
       ylim = y_range, xlab = "Year", ylab = "Precipitation (mm)", axes = TRUE)
  lines(df$year, df$predicted, col = col_pred, lty = 2)
  x_pos <- grconvertX(0.05, from = "npc", to = "user")
  y_pos <- grconvertY(0.97, from = "npc", to = "user")
  text(x_pos, y_pos, labels = title_main, adj = c(0, 1), cex = cex_title_val, font = 2)
  y_pos2 <- grconvertY(0.87, from = "npc", to = "user")
  text(x_pos, y_pos2, labels = subtitle, adj = c(0, 1), cex = cex_title_val * 0.9)
}

# Panel 1: Round A - Calibration (Early)
title_a_cal <- bquote(bold("Calibration A") ~ "(" * .(period1_years[1]) * "–" * .(period1_years[2]) * ")")
plot_panel(results_cal1, title_a_cal,
           bquote(R^2 == .(format_r2(r2_cal1))))
legend("bottomright", legend = c("Observed", "Predicted"),
       col = c("black", "red"), lty = c(1, 2), bty = "n", cex = cex_legend_val)

# Panel 2: Round A - Verification (Late)
title_a_ver <- bquote(bold("Verification A") ~ "(" * .(period2_years[1]) * "–" * .(period2_years[2]) * ")")
plot_panel(results_ver1, title_a_ver,
           bquote(R^2 == .(format_r2(r2_ver1)) ~ "|" ~ .(format_rece(re_ver1, ce_ver1))))

# Panel 3: Round B - Verification (Early)
title_b_ver <- bquote(bold("Verification B") ~ "(" * .(period1_years[1]) * "–" * .(period1_years[2]) * ")")
plot_panel(results_ver2, title_b_ver,
           bquote(R^2 == .(format_r2(r2_ver2)) ~ "|" ~ .(format_rece(re_ver2, ce_ver2))))

# Panel 4: Round B - Calibration (Late)
title_b_cal <- bquote(bold("Calibration B") ~ "(" * .(period2_years[1]) * "–" * .(period2_years[2]) * ")")
plot_panel(results_cal2, title_b_cal,
           bquote(R^2 == .(format_r2(r2_cal2))))

dev.off()

# Reset graphics parameters
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0),
    mgp = c(3, 1, 0), cex.axis = 1, cex.lab = 1, cex.main = 1.2, cex = 1)

cat(sprintf("\nDone. Plot saved to '%s'.\n", OUTPUT_PDF_FILE))
