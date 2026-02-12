# ==============================================================================
# 02_climate_modeling_utils.R
# Author: Marcos Marín-Martín
# Date: 2026-02-11
# Description:
#   Helper functions for generic climate modeling.
#   Includes robust NetCDF extraction and dynamic model fitting.
# ==============================================================================

library(ncdf4)
library(dplyr)
library(tidyr)

# ==============================================================================
# 1. NETCDF EXTRACTION FUNCTIONS
# ==============================================================================

#' Extract seasonal climate series from NetCDF (Level or Surface)
#' 
#' @param file Path to NetCDF file.
#' @param var Variable name in NetCDF.
#' @param box List with `lat` c(min, max) and `lon` c(min, max).
#' @param level Optional pressure level (e.g., 300). If NULL, assumes surface variable.
#' @param season_months Vector of months to aggregate (e.g., c(12, 1, 2) for DJF).
#' @param agg_fun Function to aggregate over space (default: mean).
#' @return Data frame with `year` and `value`.
get_netcdf_series <- function(file, var, box, level = NULL, season_months = c(12, 1, 2), agg_fun = mean) {
  
  if (!file.exists(file)) {
    warning(paste("File not found:", file))
    return(NULL)
  }
  
  nc <- tryCatch(nc_open(file), error = function(e) NULL)
  if (is.null(nc)) {
    warning(paste("Could not open NetCDF:", file))
    return(NULL)
  }
  on.exit(nc_close(nc))
  
  # --- 1. Dimensions ---
  lat <- ncvar_get(nc, "lat")
  lon <- ncvar_get(nc, "lon")
  time <- ncvar_get(nc, "time")
  
  # Smart Time Conversion (Hours vs Days)
  time_units <- ncatt_get(nc, "time", "units")$value
  if (grepl("hours", time_units, ignore.case = TRUE)) {
    dates <- as.Date(time/24, origin="1800-01-01") 
  } else {
    dates <- as.Date(time, origin="1800-01-01")
  }
  
  # --- 2. Spatial Slicing ---
  # Handle Longitude rotation if needed (0-360 vs -180-180)
  # Simple case: Assume user provides correct 0-360 coordinates
  idx_lat <- which(lat >= box$lat[1] & lat <= box$lat[2])
  idx_lon <- c(which(lon >= box$lon[1]), which(lon <= box$lon[2])) # Allows crossing 0/360 boundary if needed
  
  if (length(idx_lat) == 0 || length(idx_lon) == 0) {
    warning(paste("No grid points found in specified box for:", var))
    return(NULL)
  }
  
  # --- 3. Variable Extraction ---
  if (!is.null(level)) {
    levs <- ncvar_get(nc, "level")
    lev_idx <- which.min(abs(levs - level))
    # Check if level exists close enough
    if(abs(levs[lev_idx] - level) > 10) warning(paste("Requested level", level, "not found. Closest:", levs[lev_idx]))
    
    start <- c(1, 1, lev_idx, 1)
    count <- c(-1, -1, 1, -1)
    raw_data <- ncvar_get(nc, var, start=start, count=count)
  } else {
    raw_data <- ncvar_get(nc, var)
  }
  
  # Subset Box
  val_box <- raw_data[idx_lon, idx_lat, ]
  
  # --- 4. Spatial Aggregation ---
  # Check dimensions. If user selected 1 point, dimensions might drop.
  if (length(dim(val_box)) == 3) {
    series <- apply(val_box, 3, agg_fun, na.rm=TRUE)
  } else if (length(dim(val_box)) == 2) {
    series <- apply(val_box, 2, agg_fun, na.rm=TRUE) # Only 1 lat or lon
  } else {
    series <- val_box # Single point
  }
  
  # --- 5. Temporal Aggregation (Seasonal) ---
  df <- data.frame(
    year_raw = as.numeric(format(dates, "%Y")),
    month = as.numeric(format(dates, "%m")),
    val = series
  )
  
  # Climate Year Adjustment (e.g. Dec 2024 is Winter of 2025)
  # If season includes Dec (12) and Jan/Feb, usually we want Dec of prev year grouped with Jan of curr year.
  # We define "Climate Year" as the year of the END of the season.
  df$clim_year <- df$year_raw
  if (12 %in% season_months && any(c(1,2) %in% season_months)) {
    # If standard Winter (DJF, NDJ), shift Dec to next year
    df$clim_year <- ifelse(df$month %in% c(10, 11, 12), df$year_raw + 1, df$year_raw)
  }
  
  seasonal_series <- df %>%
    filter(month %in% season_months) %>%
    group_by(clim_year) %>%
    summarise(value = mean(val, na.rm=TRUE), n_months = n()) %>%
    filter(n_months >= length(season_months) * 0.8) %>% # Allow some missing daily days, but require months
    rename(year = clim_year)
  
  return(seasonal_series)
}


# ==============================================================================
# 2. GENERIC MODELING FUNCTIONS
# ==============================================================================

#' Fit all combinations of linear models
#' 
#' @param data Data frame containing Y (response) and X (predictors) columns.
#' @param response_var Name of the response variable column (e.g., "PCP").
#' @param predictors Vector of predictor column names.
#' @param max_vars Maximum number of variables in a single model (to prevent overfitting).
#' @return Data frame with model ranking.
fit_all_models <- function(data, response_var, predictors, max_vars = 3) {
  
  results <- data.frame(Model=character(), Formula=character(), R2_Adj=numeric(), AIC=numeric(), Num_Vars=numeric(), stringsAsFactors=FALSE)
  
  # Generate combinations
  for (k in 1:max_vars) {
    combos <- combn(predictors, k, simplify = FALSE)
    
    for (combo in combos) {
      formula_str <- paste(response_var, "~", paste(combo, collapse = " + "))
      model_name <- paste(combo, collapse = "+")
      
      tryCatch({
        mod <- lm(as.formula(formula_str), data = data)
        s <- summary(mod)
        
        results[nrow(results)+1,] <- list(
          Model = model_name,
          Formula = formula_str,
          R2_Adj = s$adj.r.squared * 100,
          AIC = AIC(mod),
          Num_Vars = k
        )
      }, error = function(e) {
        # Skip failed models (e.g. singular matrix)
      })
    }
  }
  
  # Sort by AIC (Efficiency) or R2
  results <- results %>% arrange(AIC)
  return(results)
}
