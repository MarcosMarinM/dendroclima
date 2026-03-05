# ==============================================================================
# 05_plot_calibration.R
# Author: Marcos Marín-Martín
# Date: 2026-03-05
# Description:
#   Full-Period Calibration Visualisation.
#   Generates two diagnostic plots to assess the quality of the reconstruction
#   model over the full instrumental calibration period:
#
#   1. Scatter Plot: Chronology index vs. observed precipitation, with the
#      fitted regression line, confidence band, and annotated R/R²/p statistics.
#   2. Time Series Plot: Observed vs. model-estimated precipitation overlaid,
#      allowing visual assessment of interannual tracking.
#
#   Supports two input modes:
#   - Pipeline mode: reads the _DATA.txt file from 01_main_reconstruction.R.
#   - Standalone mode: reads raw chronology + daily precipitation files and
#     computes the seasonal accumulation internally.
#
#   The input mode is selected automatically: if INPUT_FILE_RECON points to a
#   valid file, pipeline mode is used. Otherwise, standalone mode requires
#   INPUT_CHRONO_FILE and INPUT_PRECIP_FILE.
#
#   Outputs:
#   - Scatter plot (PDF + PNG).
#   - Time series plot (PDF + PNG).
# ==============================================================================

# --- 0. LIBRARIES ---
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)

# ==============================================================================
# 1. PARAMETERS
# ==============================================================================

# --- OUTPUT ---
OUTPUT_DIR <- 'PLACEHOLDER/path/to/output_dir'

# --- PIPELINE MODE: path to _DATA.txt from 01_main_reconstruction.R ---
INPUT_FILE_RECON <- 'PLACEHOLDER/path/to/reconstruction_DATA.txt'

# --- STANDALONE MODE: raw input files ---
INPUT_CHRONO_FILE <- 'PLACEHOLDER/path/to/chronology.txt'
INPUT_PRECIP_FILE <- 'PLACEHOLDER/path/to/daily_precip.txt'

# Accumulation window for standalone mode (growth year: Aug 16 to Jun 30)
WIN_START_MONTH <- 8
WIN_START_DAY   <- 16
WIN_END_MONTH   <- 6
WIN_END_DAY     <- 30

# --- MODEL OPTIONS ---
# Which chronology index to use: "res" or "std"
CHRONO_TYPE <- "res"

# Which model type for the fit: "linear" or "log"
MODEL_TYPE <- "linear"

# --- PLOT OPTIONS ---
PLOT_WIDTH  <- 8
PLOT_HEIGHT <- 5

# ==============================================================================
# 2. LOAD AND PREPARE DATA
# ==============================================================================

# Detect input mode
use_pipeline <- !grepl("PLACEHOLDER", INPUT_FILE_RECON, fixed = TRUE) &&
  file.exists(INPUT_FILE_RECON)

if (use_pipeline) {
  # --- PIPELINE MODE ---
  message("Using pipeline mode: reading ", INPUT_FILE_RECON)
  raw <- read.table(INPUT_FILE_RECON, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)

  ring_col <- CHRONO_TYPE
  rec_col  <- if (MODEL_TYPE == "linear") "rec_lineal" else "rec_log"

  required_cols <- c("year", "precipitation", ring_col, rec_col)
  missing_cols  <- setdiff(required_cols, names(raw))
  if (length(missing_cols) > 0) {
    stop("Missing columns in input file: ", paste(missing_cols, collapse = ", "))
  }

  datos_modelo <- raw |>
    select(year, precipitation, ring = all_of(ring_col),
           prec_estimada = all_of(rec_col)) |>
    filter(!is.na(precipitation), !is.na(ring))

} else {
  # --- STANDALONE MODE ---
  if (grepl("PLACEHOLDER", INPUT_CHRONO_FILE, fixed = TRUE) ||
      !file.exists(INPUT_CHRONO_FILE)) {
    stop("No valid input found. Set INPUT_FILE_RECON (pipeline) or ",
         "INPUT_CHRONO_FILE + INPUT_PRECIP_FILE (standalone).")
  }
  message("Using standalone mode: reading raw chronology + precipitation files.")

  crono <- read.table(INPUT_CHRONO_FILE, header = TRUE)
  prec  <- read.table(INPUT_PRECIP_FILE, header = TRUE)

  # Process daily precipitation into growth-year accumulation
  prec$date <- ymd(prec$date)

  prec <- prec |>
    mutate(
      year  = year(date),
      month = month(date),
      day   = day(date),
      growth_year = case_when(
        month == WIN_START_MONTH & day >= WIN_START_DAY ~ year + 1L,
        month > WIN_START_MONTH                         ~ year + 1L,
        month < WIN_END_MONTH                           ~ year,
        month == WIN_END_MONTH & day <= WIN_END_DAY     ~ year,
        TRUE                                            ~ NA_integer_
      )
    )

  prec_acum <- prec |>
    filter(!is.na(growth_year)) |>
    group_by(growth_year) |>
    summarise(precipitation = sum(precipitation, na.rm = TRUE), .groups = "drop")

  datos_modelo <- crono |>
    inner_join(prec_acum, by = c("year" = "growth_year")) |>
    filter(!is.na(.data[[CHRONO_TYPE]])) |>
    rename(ring = all_of(CHRONO_TYPE))

  # Fit model to get estimated values
  if (MODEL_TYPE == "linear") {
    mod <- lm(precipitation ~ ring, data = datos_modelo)
    datos_modelo$prec_estimada <- predict(mod)
  } else {
    mod <- lm(log(precipitation) ~ ring, data = datos_modelo)
    datos_modelo$prec_estimada <- exp(predict(mod))
  }
}

if (nrow(datos_modelo) < 10) {
  stop("Less than 10 overlapping years with observed precipitation found.")
}

# ==============================================================================
# 3. COMPUTE STATISTICS
# ==============================================================================

if (MODEL_TYPE == "linear") {
  modelo <- lm(precipitation ~ ring, data = datos_modelo)
} else {
  modelo <- lm(log(precipitation) ~ ring, data = datos_modelo)
}

r_val  <- cor(datos_modelo$precipitation, datos_modelo$prec_estimada,
              use = "complete.obs")
R2_val <- summary(modelo)$r.squared
p_val  <- summary(modelo)$coefficients[2, 4]

p_text     <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
stats_text <- sprintf("r = %.2f | r\u00b2 = %.2f | %s", r_val, R2_val, p_text)

# Create output directory
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ==============================================================================
# 4. SCATTER PLOT
# ==============================================================================

plot_scatter <- ggplot(datos_modelo, aes(x = ring, y = precipitation)) +
  geom_point(color = "black", size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", color = "red", linetype = "dashed",
              se = TRUE, fill = "gray80") +
  labs(
    title    = paste0("Calibration: ", toupper(CHRONO_TYPE),
                      " vs Precipitation (", MODEL_TYPE, " model)"),
    subtitle = stats_text,
    x        = "Chronology index",
    y        = "Precipitation (mm)"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

ggsave(file.path(OUTPUT_DIR, "calibration_scatter.pdf"),
       plot = plot_scatter, width = PLOT_WIDTH, height = PLOT_HEIGHT)
ggsave(file.path(OUTPUT_DIR, "calibration_scatter.png"),
       plot = plot_scatter, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = 300)

# ==============================================================================
# 5. TIME SERIES PLOT
# ==============================================================================

datos_long <- datos_modelo |>
  select(year, Observed = precipitation, Estimated = prec_estimada) |>
  pivot_longer(cols = c(Observed, Estimated),
               names_to = "Serie", values_to = "Precipitation") |>
  mutate(Serie = factor(Serie, levels = c("Observed", "Estimated")))

plot_timeseries <- ggplot(datos_long,
                          aes(x = year, y = Precipitation,
                              color = Serie, linetype = Serie)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(name = "",
                     values = c("Observed" = "black", "Estimated" = "red")) +
  scale_linetype_manual(name = "",
                        values = c("Observed" = "solid", "Estimated" = "dashed")) +
  labs(
    title    = "Calibration Period: Observed vs Estimated",
    subtitle = stats_text,
    x        = "Year",
    y        = "Precipitation (mm)"
  ) +
  theme_classic() +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle    = element_text(hjust = 0.5, size = 12),
    legend.position  = c(0.85, 0.88),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.width = unit(1.5, "cm")
  )

ggsave(file.path(OUTPUT_DIR, "calibration_timeseries.pdf"),
       plot = plot_timeseries, width = PLOT_WIDTH, height = PLOT_HEIGHT)
ggsave(file.path(OUTPUT_DIR, "calibration_timeseries.png"),
       plot = plot_timeseries, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = 300)

cat(sprintf("Plots saved to: %s\n", OUTPUT_DIR))

