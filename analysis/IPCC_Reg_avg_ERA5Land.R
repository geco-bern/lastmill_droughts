########### Script for data post-processing to create regional averages (see IPCC regions, AR6) ###########
### load relevant libraries:
library(dplyr)
library(rJava)
library(loadeR.java)
library(transformeR)
library(loadeR)
library(visualizeR)
library(geoprocessoR)
library(terra)
library(ncdf4)
library(rnaturalearth)
library(sf)
library(here)
library(sp)

# Define what to run:
flag_only_run_CHregion         <- TRUE  # if TRUE, only runs CH, otherwise all regions (including CH)

#flag_only_run_2025_2026_update <- FALSE # if TRUE, only runs 2025/2026, otherwise only -2024
flag_only_run_2025_2026_update <- TRUE # if TRUE, only runs 2025/2026, otherwise only -2024


# 1. Define FILEPATH as a string (Climate4R needs a string path)
if (flag_only_run_2025_2026_update){
  nc_path <- here("data/new_ERA5_PCWD_ANNMAX.nc")      # for 2025,2026
} else {
  nc_path <- here("data/ERA5Land_cons_PCWD_ANNMAX.nc") # for paper
}

## read metadata variables using a temporary connection
nc_pwcd  <- nc_open(nc_path)
lon      <- ncvar_get(nc_pwcd, varid="lon")
lat      <- ncvar_get(nc_pwcd, varid="lat")
time     <- ncvar_get(nc_pwcd, varid="time")

reference_date  <- as.Date("2001-01-01")
time_dates_1850 <- reference_date + time
nc_close(nc_pwcd)


# 2. Land-Sea Mask Setup --------------------------------------------------
land <- ne_countries(scale = "medium", returnclass = "sf")
land <- st_transform(land, crs = "+proj=longlat +datum=WGS84")

grid_df <- expand.grid(lon = lon, lat = lat)
grid_sf <- st_as_sf(grid_df, coords = c("lon", "lat"), crs = st_crs(land))

grid_sf$in_land <- st_intersects(grid_sf, land, sparse = FALSE) %>% rowSums() > 0

# Match matrix orientation
land_sea_mask <- matrix(grid_sf$in_land, nrow = length(lon), ncol = length(lat), byrow = FALSE)
land_sea_mask <- t(land_sea_mask)


# 3. Load Climate Data ONCE (Pull out of the loop for speed!) -----------
grid <- loadGridData(dataset = nc_path, var = "pcwd_annmax")

# Apply land mask directly to the full grid object
data_array   <- grid$Data
masked_array <- data_array

for (t in 1:dim(data_array)[1]) {
  # Apply mask to each time slice matrix
  slice <- masked_array[t, , ]
  slice[!land_sea_mask] <- NA
  masked_array[t, , ] <- slice
}
grid$Data <- masked_array


# 4. Reference Regions Setup ---------------------------------------------
source(here::here("analysis","define_reference_regions.R"))
refregion_objects <- define_reference_regions()

regions <- names(refregion_objects)
# refregion_objects['WCE']
# refregion_objects['WCE_W'] # new regions
# refregion_objects['WCE_W'] # new regions
# refregion_objects['CHE'] # new regions


grid <- setGridProj(grid = grid, proj = proj4string(refregion_objects[[1]]))


########## for ERA5Land ###########################
regional_results <- list()

# Get the total number of time steps from your main grid object
total_time_steps <- dim(grid$Data)[1]

if (flag_only_run_CHregion){
  regions_to_loop <- c("CHE","WCE_W","WCE") # only subset
} else {
  regions_to_loop <- regions # all, including CHE, WCE_W, etc..
}
for (region in regions_to_loop) {
  message("Processing region: ", region)

  region_object <- refregion_objects[[region]]
  if (is.null(region_object)) {
    warning(paste("Region not found:", region))
    next
  }

  # Crop grid to specific IPCC region using Climate4R overlay
  grid_region <- tryCatch({
    overGrid(grid, region_object)
  }, error = function(e) {
    NULL
  })

  # CRITICAL SAFETY CHECK: If overGrid failed or returned nothing
  if (is.null(grid_region) || is.null(grid_region$Data)) {
    warning(paste("No overlapping grid points found for region:", region))
    regional_results[[region]] <- rep(NA_real_, total_time_steps)
    next
  }

  region_data <- grid_region$Data
  data_dims   <- dim(region_data)

  # SAFETY CHECK: ensure region is large enough:
  n_valid_cells <- max(
    apply(!is.na(region_data), 1, sum),
    na.rm = TRUE
  )

  message(region, ": ", n_valid_cells, " contributing grid cells")

  if (!is.finite(n_valid_cells) || n_valid_cells == 0) {
    regional_results[[region]] <- rep(NA_real_, total_time_steps)
    next
  }

  # CRITICAL SAFETY CHECK 2: Ensure the data is actually 3D [Time, Lat, Lon]
  # If it collapsed to less than 3 dimensions, there are no valid grid cells.
  if (length(data_dims) < 3 || any(data_dims[2:3] == 0)) {
    warning(paste("Spatial dimensions collapsed for region:", region, "- filling with NA"))
    regional_results[[region]] <- rep(NA_real_, total_time_steps)
    next
  }

  # Compute regional spatial average across spatial dimensions (dims 2 and 3)
  spatial_avg <- apply(region_data, 1, function(slice) {
    # If the entire slice is NA, mean(..., na.rm=TRUE) returns NaN.
    # Let's convert NaN to a clean NA_real_
    res <- mean(slice, na.rm = TRUE)
    if (is.nan(res)) return(NA_real_)
    return(res)
  })

  # Store the direct numeric vector into your tracking list
  regional_results[[region]] <- spatial_avg
}
# OUTPUT WAS:
#     [...]
#     SAU: 45 contributing grid cells
#     Processing region: NZ
#     NZ: 10 contributing grid cells
#     Processing region: WCE_E
#     WCE_E: 93 contributing grid cells
#     Processing region: WCE_W
#     WCE_W: 31 contributing grid cells
#     Processing region: CHE
#     CHE: 1 contributing grid cells
#     Processing region: BERN (25km)
#     BERN (25km): 0 contributing grid cells
#     Processing region: Bern/Zollikofen
#     Bern/Zollikofen: 0 contributing grid cells

#save calculated list (note run once for base data and once again for 25_26 data)
if (flag_only_run_CHregion){
  if (flag_only_run_2025_2026_update){
    saveRDS(regional_results, file=(here("data/regionalResults_CH_only_ERA5Land_cons_25_26.RData"))) # for 2025,2026,newRegion
  } else {
    saveRDS(regional_results, file=(here("data/regionalResults_CH_only_ERA5Land_cons.RData")))       # for old-period,newRegion
  }
} else {
  if (flag_only_run_2025_2026_update){
    saveRDS(regional_results, file=(here("data/regionalResults_ERA5Land_cons_25_26.RData"))) # for 2025,2026
  } else {
    saveRDS(regional_results, file=(here("data/regionalResults_ERA5Land_cons.RData")))       # for paper
  }
}



## check grid alignment
# Load relevant libraries
library(dplyr)
library(rJava)
library(loadeR.java)
library(transformeR)
library(loadeR)
library(visualizeR)
#library(geoprocessoR)
library(terra)
library(ncdf4)
library(rnaturalearth)
library(sf)
library(raster)
library(ggplot2)
library(viridis)
library(here)

# 1. Define FILEPATH as a string (Fixes loadGridData crash)
nc_path <- here("data/ERA5Land_cons_PCWD_ANNMAX.nc")

## Read metadata spatial variables
nc_pwcd <- nc_open(nc_path)
lon     <- ncvar_get(nc_pwcd, varid="lon")
lat     <- ncvar_get(nc_pwcd, varid="lat")
nc_close(nc_pwcd)

# 2. Build the Land Mask Grid Matrix
land    <- ne_countries(scale = "medium", returnclass = "sf")
land    <- st_transform(land, crs = "+proj=longlat +datum=WGS84")
grid_df <- expand.grid(lon = lon, lat = lat)
grid_sf <- st_as_sf(grid_df, coords = c("lon", "lat"), crs = st_crs(land))

grid_sf$in_land <- st_intersects(grid_sf, land, sparse = FALSE) %>% rowSums() > 0

# Match matrix orientation: Row = Lon, Col = Lat
land_sea_mask <- matrix(grid_sf$in_land, nrow = length(lon), ncol = length(lat), byrow = FALSE)
land_sea_mask <- t(land_sea_mask) # Row = Lat, Col = Lon to match climate4R structure

# 3. Load Data Using Climate4R (Pass file path string)
grid_test       <- loadGridData(dataset = nc_path, var = "pcwd_annmax")
grid_array_test <- grid_test$Data

# 4. Pick step 2 slice and apply the mask directly
# climate4R matrix dimensions: [Time, Lat, Lon]
plot_data <- grid_array_test[2, , ]
plot_data[!land_sea_mask] <- NA # Set ocean pixels to NA

# 1. Try passing the plot_data WITHOUT transposing it first
# (Or if it was already transposed, try the reverse of what you had)
r <- raster(plot_data, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
crs(r) <- "+proj=longlat +datum=WGS84"

# 2. Check the flip direction: NetCDF matrices are usually inverted vertically
r <- flip(r, direction='y')

# 6. Prepare Dataframe for ggplot
r_df <- as.data.frame(r, xy = TRUE)
colnames(r_df) <- c("x", "y", "layer")

# 7. Render Spatial Map Verification
ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = layer), na.rm = TRUE) +
  scale_fill_viridis_c(option = "C", na.value = "transparent") +
  labs(title = "Annmax PCWD in 1952 (Spatial Check)", x = "Longitude", y = "Latitude") +
  theme_classic() +
  geom_sf(data = land, fill = NA, color = "black", lwd = 0.4) +
  coord_sf(xlim = range(lon), ylim = range(lat), expand = FALSE) +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  guides(fill = guide_colorbar(title = "PCWD (mm)"))


# ### Compare bilinear and conservatively interpolated PCWD
#
# regional_results_ERA5Land_bil <- readRDS(here("data/regionalResults_ERA5Land.RData")) #1950 - 2024   # TODO: missing data/regionalResults_ERA5Land.RData
# ### to exclude potential precipitation problems in ERA5-Land only include years from 1970s onwards
# regional_results_ERA5Land_bil <- lapply(regional_results_ERA5Land_bil, function(x) x[21:75])
#
#
# regional_results_ERA5Land_con <- readRDS(here("data/regionalResults_ERA5Land_cons.RData")) #1950 - 2024 # TODO: create data/regionalResults_ERA5Land_cons.RData
# ### to exclude potential precipitation problems in ERA5-Land only include years from 1970s onwards
# regional_results_ERA5Land_con <- lapply(regional_results_ERA5Land_con, function(x) x[21:75])
#
#
# library(ggplot2)
# library(tidyr)
# library(dplyr)
#
# # 1. Define the correct year sequence (indices 21:75 map to 1970:2024)
# years_vector <- 1970:2024
#
# # 2. Extract the MED vectors from both lists
# med_bil <- regional_results_ERA5Land_bil[["CAU"]]
# med_con <- regional_results_ERA5Land_con[["CAU"]]
#
# # 3. Create a combined dataframe structured for ggplot
# plot_df <- data.frame(
#   Year = years_vector,
#   Bilinear = med_bil,
#   Conservative = med_con
# ) |>
#   # Reshape data into long format for clean plotting aesthetics
#   pivot_longer(
#     cols = c(Bilinear, Conservative),
#     names_to = "Remapping_Method",
#     values_to = "PCWD"
#   )
#
# # 4. Generate the Time-Series Comparison Plot
# ggplot(plot_df, aes(x = Year, y = PCWD, color = Remapping_Method, linetype = Remapping_Method)) +
#   geom_line(lwd = 1) +
#   geom_point(size = 1.5) +
#   # Custom aesthetic layout adjustments
#   scale_color_manual(values = c("Bilinear" = "#1f77b4", "Conservative" = "#ff7f0e")) +
#   scale_linetype_manual(values = c("Bilinear" = "solid", "Conservative" = "dashed")) +
#   labs(
#     title = "ERA5-Land PCWD Comparison: Central Australia Region (CAU)",
#     subtitle = "Comparing Bilinear vs. Conservative Remapping Schemes (1970–2024)",
#     x = "Year",
#     y = "Potential Cumulative Water Deficit (mm)",
#     color = "Remapping Method",
#     linetype = "Remapping Method"
#   ) +
#   theme_classic(base_size = 13) +
#   theme(
#     legend.position = "bottom",
#     plot.title = element_text(face = "bold", hjust = 0.5),
#     plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
#     panel.grid.major.y = element_line(color = "gray90") # Easier to read historical variance
#   ) +
#   scale_x_continuous(breaks = seq(1970, 2024, by = 5))
