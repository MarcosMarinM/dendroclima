# dendroclima: dendroclimatic analysis pipeline

This repository contains a structured R project for processing, analysing, and visualising dendroclimatic data. The workflow is organised into six main stages, designed to guide the user from raw data to advanced climate modeling.

## Project structure

The project is now organised into the following directories, reflecting a modular analytical workflow:

*   **`00_utils`**:
    *   Contains helper scripts and shared functions used across the project (e.g., `dendroTools`).
    *   **Key scripts**: `00_utils_dendro.R`, `01_dendro_tools.R`.

*   **`01_dendro_processing`**:
    *   **Goal**: Process raw tree-ring width (RWL) files into standardised chronologies.
    *   **Key scripts**: `02_main_detrending_analysis.R`
    *   **Input**: Raw `.rwl` files.
    *   **Output**: Standardised (`.std`) and Residual (`.res`) chronologies, Quality Control Plots.

*   **`02_climate_extraction`**:
    *   **Goal**: Extract raw climate data from external sources (TerraClimate, PHYDA, SLP).
    *   **Key scripts**: `01_extract_phyda.R`, `04_extract_terra_point.R`.

*   **`03_climate_preparation`**:
    *   **Goal**: Format separate text files and calculate derived climate indices.
    *   **Key scripts**: `01_format_rocio_txt.R`, `03_calc_spei.R`, `08_sweep_seasonal_nao.R`.
    *   **Output**: Cleaned `.txt` files ready for correlation analysis.

*   **`04_climate_analysis`**:
    *   **Goal**: Analyze relationships between tree-ring growth and climate variables.
    *   **Key scripts**: 
        *   **`00_dendro_hybrid_analysis.R`**: Standalone script for quick correlation analysis, cross-validation, and climate window detection.
        *   `01_main_correlations.R`: Core engine for batch correlation of all chronologies.
        *   `02_loop_rwl_clim.R`, `03_loop_rwl_clim_parallel.R`: Standardization and correlation loops.
    *   **Input**: Chronologies from `01` + climate data from `03`.
    *   **Output**: Correlation heatmaps, statistical reports.
    *   **Python tools**: Includes scripts (`08_spatial_corr.py` to `12_spatial_corr_eobs_crop.py`) for spatial correlation analysis.

*   **`05_reconstruction`**:
    *   **Goal**: Reconstruct past climate variability based on tree-ring/climate relationships.
    *   **Key scripts**: `01_main_reconstruction.R`
    *   **Input**: Calibrated chronologies, instrumental climate Records.
    *   **Output**: Reconstructed Time Series (Text & Plots), Verification Statistics (R2, RE, CE).

*   **`06_climate_modeling`**:
    *   **Goal**: Explore drivers of precipitation using atmospheric circulation models.
    *   **Key scripts**: `01_modeling_precip_battle_royale.R`
    *   **Input**: Reanalysis Data (NCEP/NCAR).
    *   **Output**: Model Rankings (AIC/R2), Physical Mechanism Analysis.

## Prerequisites

*   **R** (>= 4.0.0)
*   **RStudio** (recommended)
*   **Python** (for spatial correlations) with `numpy`, `netCDF4`, `matplotlib`, `cartopy`, `scipy`.

### R packages

Ensure the following packages are installed:

```r
install.packages(c(
  "tidyverse", "dplR", "dendroTools", "SPEI", "ggplot2", 
  "reshape2", "future", "future.apply", "ncdf4", "sf", 
  "terra", "patchwork", "zoo", "colourvalues", "magick"
))
```

## Usage

1.  **Setup**: Open a project in RStudio.
6.  **Configuration**: 
    *   Navigate to the script you wish to run.
    *   Locate the **PARAMETERS** section at the top of the script.
    *   **IMPORTANT**: All file paths have been replaced with `PLACEHOLDER/path/to/...`. You **MUST** update these paths to match your local system structure before running any script.
    *   Libraries are now consolidated at the top of each script. Ensure you have all required packages installed.
7.  **Execution**: Run scripts sequentially following the directory numbers (01 -> 06).