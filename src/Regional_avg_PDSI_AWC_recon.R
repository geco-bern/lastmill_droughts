#!/usr/bin/env Rscript
# =============================================================================
# Script:  run_PDSI_regional_averaging.R
# Author:  Patricia Helpap
# Date:    Sys.Date()
# Purpose: Compute IPCC AR6 regional averages of scPDSI from AWC reconstruction
#          (ModE-Sim ensemble, 1420-2009). Designed to be run in a terminal /
#          screen session:
#
#   screen -S pdsi
#   Rscript run_PDSI_regional_averaging.R > logs/pdsi_regional.log 2>&1
#
# =============================================================================

# ── 0. Libraries ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(slider)
  library(lubridate)
  library(patchwork)
  library(ggplot2)
  library(ncdf4)
  library(reshape2)
  library(ggpubr)
  library(ggtext)
  library(maps)
  library(terra)
  library(rnaturalearth)
  library(sf)
  library(sp)
  library(RColorBrewer)
  library(tidyverse)
  library(abind)
  library(gridExtra)
  library(ggsci)
  library(purrr)
  library(cowplot)
  library(extRemes)
  library(rJava)
  library(loadeR.java)
  library(transformeR)
  library(loadeR)
  library(visualizeR)
})

cat("Libraries loaded.\n")

# ── 1. Read lon/lat from one test file ───────────────────────────────────────
PDSI_recon_01 <- "/Users/phelpap/Documents/lastmilldroughts_data/scPDSI_AWC_recon/ModE-Sim_m001_PDSI_1420-2009.nc"

nc_test        <- nc_open(PDSI_recon_01)
var            <- ncvar_get(nc_test, varid = "pdsi")
lon            <- ncvar_get(nc_test, varid = "longitude")
lat            <- ncvar_get(nc_test, varid = "latitude")
time           <- ncvar_get(nc_test, varid = "time")
nc_close(nc_test)

# Convert time (days since 1970-01-01) to dates
reference_date <- as.Date("1970-01-01")
time_dates     <- reference_date + time
cat(sprintf("Time axis: %s  →  %s  (%d time steps)\n",
            format(min(time_dates)), format(max(time_dates)), length(time_dates)))

# ── 2. Load IPCC AR6 reference regions ───────────────────────────────────────
load(here::here("data", "Reference_regions", "IPCC-WGI-reference-regions-v4_R.rda"))
refregions <- as(IPCC_WGI_reference_regions_v4, "SpatialPolygons")

cat("IPCC reference regions loaded.\n")

# ── 3. Configuration ─────────────────────────────────────────────────────────

# All land regions (uncomment to run the full set):
regions <- c("GIC", "NWN", "NEN", "WNA", "CNA", "ENA", "NCA", "SCA", "CAR",
             "NWS", "NSA", "NES", "SAM", "SWS", "SES", "SSA", "NEU", "WCE",
             "EEU", "MED", "SAH", "WAF", "CAF", "NEAF", "SEAF", "WSAF", "ESAF",
             "MDG", "RAR", "WSB", "ESB", "RFE", "WCA", "ECA", "TIB", "EAS",
             "ARP", "SAS", "SEA", "NAU", "CAU", "EAU", "SAU", "NZ")
#regions <- c("GIC", "MED", "WNA", "ARP", "CAU")   # testing — remove / expand as needed

path  <- "/Users/phelpap/Documents/lastmilldroughts_data/scPDSI_AWC_recon"
files <- list.files(path, full.names = TRUE, pattern = "m[0-9]{3}")

# Derive ensemble member IDs directly from filenames (avoids the NA-column bug)
ensemble_members <- regmatches(basename(files), regexpr("m[0-9]{3}", basename(files)))

cat(sprintf("Found %d file(s): %s\n", length(files), paste(ensemble_members, collapse = ", ")))

# ── 4. Main loop: regions × ensemble members ──────────────────────────────────
regional_results_AWC_recon <- list()

for (region in regions) {

  cat(sprintf("\n── Region: %s ──\n", region))

  region_object <- refregions[c(region)]

  if (is.null(region_object)) {
    warning(paste("Region not found:", region))
    next
  }

  bbox           <- bbox(region_object)
  spatial_avg_list <- list()

  for (i in seq_along(files)) {
    file <- files[i]
    em   <- ensemble_members[i]

    cat(sprintf("  Processing %s ...\n", em))

    if (!file.exists(file)) {
      warning(paste("Missing file for ensemble member:", em))
      next
    }

    # Load grid data
    grid <- loadGridData(dataset = file, var = "pdsi")

    # Set spatial projection to match reference regions
    grid <- setGridProj(grid = grid, proj = proj4string(refregions))

    # Crop to bounding box of the region
    grid_region <- subsetGrid(grid,
                              lonLim  = c(bbox["x", "min"], bbox["x", "max"]),
                              latLim  = c(bbox["y", "min"], bbox["y", "max"]),
                              outside = FALSE)

    # Build a mask for the exact polygon shape (not just the bounding box)
    coords    <- getCoordinates(grid_region)
    points_sp <- SpatialPoints(expand.grid(coords$x, coords$y),
                               proj4string = CRS(proj4string(refregions)))
    inside         <- !is.na(over(points_sp, region_object))
    inside_matrix  <- matrix(inside,
                             nrow = length(coords$x),
                             ncol = length(coords$y))
    inside_matrix  <- t(inside_matrix)   # → [lat × lon] matches Data[t, lat, lon]

    for (t in 1:dim(grid_region$Data)[1]) {
      grid_region$Data[t, , ][!inside_matrix] <- NA
    }

    # Spatial average ignoring masked NAs
    spatial_avg <- apply(grid_region$Data, 1, function(slice) {
      mean(slice, na.rm = TRUE)
    })
    spatial_avg_list[[em]] <- spatial_avg  # keyed by "m001", "m002", etc.

    cat(sprintf("    → %d time steps, mean PDSI = %.3f\n",
                length(spatial_avg), mean(spatial_avg, na.rm = TRUE)))
  }

  # Combine ensemble members as columns
  result_array_region                    <- do.call(cbind, spatial_avg_list)
  regional_results_AWC_recon[[region]]   <- result_array_region
}

cat("\nAll regions processed.\n")

# ── 5. Save output ────────────────────────────────────────────────────────────
out_path <- here::here("data", "regionalResults_scPDSI_1420_2009_AWCrecon_all.RData")
saveRDS(regional_results_AWC_recon, file = out_path)
cat(sprintf("Results saved to: %s\n", out_path))

# ── 6. Quick diagnostic plot for GIC / m001 ───────────────────────────────────
if ("MED" %in% names(regional_results_AWC_recon) &&
    "m001" %in% colnames(regional_results_AWC_recon[["GIC"]])) {

  pdsi_ts   <- regional_results_AWC_recon[["GIC"]][, "m001"]
  n_steps   <- length(pdsi_ts)
  time_axis <- seq(1420, by = 1/12, length.out = n_steps)

  df_plot <- data.frame(time = time_axis, pdsi = pdsi_ts)

  p <- ggplot(df_plot, aes(x = time, y = pdsi)) +
    geom_line(colour = "steelblue", linewidth = 0.4) +
    geom_smooth(method = "loess", span = 0.1,
                colour = "red", linewidth = 0.8, se = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    labs(title = "GIC Regional PDSI — m001 (AWC reconstruction)",
         x = "Year", y = "PDSI") +
    theme_classic()

  plot_path <- here::here("data", "MED_m001_PDSI_timeseries.png")
  ggsave(plot_path, plot = p, width = 12, height = 4, dpi = 150)
  cat(sprintf("Diagnostic plot saved to: %s\n", plot_path))
}

cat("\nDone.\n")
