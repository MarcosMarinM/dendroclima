# ==============================================================================
# 00_utils_dendro.R
# Author: Marcos Marín-Martín
# Date: 2026-02-09
# Description:
#   General Utility Functions for Dendroclimatology.
#   A collection of helper functions used across the project for common tasks 
#   such as data formatting, simple statistical calculations, and file handling.
#
#   Contents:
#   - clean_dates(): Standardizes date formats.
#   - load_rwl_safe(): Wrapper for reading RWL files with error handling.
#   - basic_stats(): Computes mean, sd, se, etc., with NA handling.
#   - ... and other miscellaneous helpers.
#
#   Usage:
#   - Source this file at the beginning of analysis scripts to access these tools.
# ==============================================================================


library(lubridate)
library(dplyr)
library(tidyr)
library(stringr)

# --- Climate Data Processing ---

#' Read daily climate data from TXT file
#' 
#' @param file_path Path to the climate data file
#' @param var_name Name of the variable to extract/rename
#' @return A dataframe with year, month, day, date, and value
read_daily_climate_data <- function(file_path, var_name) {
  cat("Reading:", file_path, "\n")
  first_line <- readLines(file_path, n = 1)
  sep <- ifelse(grepl("\t", first_line), "\t", ifelse(grepl(" ", first_line), "", NA))
  if (is.na(sep)) stop(paste("Could not determine separator for:", file_path))
  
  data <- tryCatch(read.table(file_path, header = TRUE, sep = sep, stringsAsFactors = FALSE, na.strings=c("NA", "-9999", "9999")),
                   error = function(e) stop(paste("Error reading", file_path, ":", e$message)))
  
  if (!("date" %in% colnames(data))) stop(paste("'date' column not found in", file_path))
  if (!(var_name %in% colnames(data))) {
    if (ncol(data) == 2) {
      actual_var_col <- setdiff(colnames(data), "date")
      warning(paste0("'", var_name, "' not found in ", file_path, ". Using '", actual_var_col, "' instead."))
      colnames(data)[colnames(data) == actual_var_col] <- var_name
    } else stop(paste("'", var_name, "' not found and >2 columns in", file_path))
  }
  
  data$date_obj <- ymd(as.character(data$date))
  if(any(is.na(data$date_obj))) {
    warning(paste(sum(is.na(data$date_obj)), "dates not parsed in", basename(file_path)))
  }
  data %>% filter(!is.na(date_obj)) %>%
    mutate(year = year(date_obj), month = month(date_obj), day = day(date_obj)) %>%
    select(year, month, day, date = date_obj, value = !!sym(var_name))
}

#' Aggregate daily climate data to monthly, seasonal, and annual
#' 
#' @param daily_data Dataframe returned by read_daily_climate_data
#' @param agg_fun Aggregation function ("mean" or "sum")
#' @return List of dataframes (monthly, seasonal, annual)
aggregate_climate_data <- function(daily_data, agg_fun = "mean") {
  if (nrow(daily_data) == 0) return(list(monthly=data.frame(), seasonal=data.frame(), annual=data.frame()))
  
  fn <- if (agg_fun == "sum") sum else mean
  
  monthly_agg <- daily_data %>% filter(!is.na(value)) %>% group_by(year, month) %>%
    summarise(value = fn(value, na.rm = TRUE), n_days = n(), .groups = "drop")
  
  annual_agg <- daily_data %>% filter(!is.na(value)) %>% group_by(year) %>%
    summarise(value = fn(value, na.rm = TRUE), n_days = n(), .groups = "drop")
  
  seasonal_data_prep <- daily_data %>% filter(!is.na(value)) %>%
    mutate(season_year = ifelse(month == 12, year + 1, year),
           season = case_when(month %in% c(12, 1, 2) ~ "DJF", month %in% c(3, 4, 5) ~ "MAM",
                              month %in% c(6, 7, 8) ~ "JJA", month %in% c(9, 10, 11) ~ "SON",
                              TRUE ~ NA_character_)) %>%
    filter(!is.na(season))
  
  seasonal_agg <- seasonal_data_prep %>% group_by(season_year, season) %>%
    summarise(value = fn(value, na.rm = TRUE), n_days = n(), .groups = "drop") %>%
    rename(year = season_year)
  
  list(monthly = select(monthly_agg, -any_of("n_days")),
       seasonal = select(seasonal_agg, -any_of("n_days")),
       annual = select(annual_agg, -any_of("n_days")))
}

# --- Plotting Helpers ---

#' Calculate dynamic Y-axis limits and breaks for plotting
#' 
#' @param correlations Vector of correlation values
#' @return List with 'limits' and 'breaks'
get_dynamic_y_limits_and_breaks <- function(correlations) {
  valid_correlations <- correlations[!is.na(correlations) & is.finite(correlations)]
  if (length(valid_correlations) == 0) return(list(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)))
  
  max_abs_corr <- max(abs(valid_correlations), na.rm = TRUE)
  if (max_abs_corr == 0) return(list(limits = c(-0.1, 0.1), breaks = seq(-0.1, 0.1, 0.05)))
  
  if (max_abs_corr <= 0.1) { upper_limit <- 0.1; interval <- 0.05 }
  else if (max_abs_corr <= 0.25) { upper_limit <- 0.25; interval <- 0.05 }
  else if (max_abs_corr <= 0.5) { upper_limit <- 0.5; interval <- 0.1 }
  else if (max_abs_corr <= 0.75) { upper_limit <- 0.75; interval <- 0.25 }
  else { upper_limit <- 1.0; interval <- 0.25 }
  
  upper_limit <- min(upper_limit, 1.0)
  limits <- c(-upper_limit, upper_limit)
  breaks <- unique(round(seq(limits[1], limits[2], by = interval), 2))
  
  return(list(limits = limits, breaks = breaks))
}

#' Parse information from chronology filename
#' 
#' @param filename_txt Filename string (e.g. "site_std_mid.txt")
#' @return List with site, type, age_group, and original_filename
get_info_from_chron_filename <- function(filename_txt) {
  base_name <- sub("\\.txt$", "", filename_txt)
  parts <- unlist(strsplit(base_name, "_"))
  type <- "UNKNOWN"; age_group <- "UNKNOWN"; site <- if(length(parts)>=1) parts[1] else "UNKNOWN"
  
  if (length(parts) > 1) {
    potential_age <- tolower(parts[length(parts)])
    if (potential_age %in% c("old", "mid", "yng")) {
      age_group <- potential_age
      type <- if (length(parts) > 2) paste(parts[2:(length(parts)-1)], collapse="_") else site
    } else {
      type <- paste(parts[2:length(parts)], collapse="_"); age_group <- "all"
    }
  } else { type <- base_name; age_group <- "all" }
  
  type_norm <- toupper(type)
  if (type_norm %in% c("ΔBI", "DBI", "DBI")) type_norm = "DBI" 
  return(list(site=site, type=type_norm, age_group=age_group, original_filename=filename_txt))
}
