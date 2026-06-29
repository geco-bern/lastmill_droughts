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

#lat lon and time info from netcdf files
##read in data
input_file_1850 <- "/storage/research/giub_geco/data_2/scratch/phelpap/ModESim/m001_tidy/04_result_1850/PCWD_ANNMAX.nc"
input_file_1420 <- "/storage/research/giub_geco/data_2/scratch/phelpap/ModESim/m001_tidy/04_result_1420/PCWD_ANNMAX.nc"

nc_pwcd_1850 <- nc_open(input_file_1850)
pcwd_annmax_1850 = ncvar_get(nc_pwcd_1850, varid="pcwd_annmax")
lon = ncvar_get(nc_pwcd_1850, varid="lon")
lat = ncvar_get(nc_pwcd_1850, varid="lat")
time_1850 = ncvar_get(nc_pwcd_1850, varid="time")
# Convert to actual dates (days since 2001-01-01)
reference_date <- as.Date("2001-01-01")
time_dates_1850 <- reference_date + time_1850

# # Print the resulting dates
# print(time_dates)

nc_close(nc_pwcd_1850)

#1420 file for 1420 time
nc_pwcd_1420 <- nc_open(input_file_1420)
pcwd_annmax_1420 = ncvar_get(nc_pwcd_1420, varid="pcwd_annmax")
time_1420 = ncvar_get(nc_pwcd_1420, varid="time")
# Convert to actual dates (days since 2001-01-01)
time_dates_1420 <- reference_date + time_1420

# # Print the resulting dates
# print(time_dates)

nc_close(nc_pwcd_1420)


# ### read in ModE-Sim land-sea mask information (binary):
# ### land mask included?
# input_file <- "/storage/research/giub_geco/data_2/scratch/phelpap/ModESim/ModESim_forcings/stable_inputdata/T63GR15_jan_surf.nc"
# nc <- nc_open(input_file)
# slm = ncvar_get(nc, varid="SLM") # 1 = land, 0 = ocean
# nc_close(nc)
# slm <- t(slm)  # Transpose the land-sea mask to match dimensions of data


# Step 1: Create a land-sea mask from the natural earth dataset
land <- ne_countries(scale = "medium", returnclass = "sf")
land <- st_transform(land, crs = "+proj=longlat +datum=WGS84")

# Step 2: Create a grid with the same lon/lat as gini_values_1850
grid_df <- expand.grid(lon = lon, lat = lat)

# Step 3: Check which grid cells are land or ocean based on the land polygons
grid_sf <- st_as_sf(grid_df, coords = c("lon", "lat"), crs = st_crs(land))

# Step 4: Use `st_intersects` to check if each grid cell is land (1) or ocean (0)
grid_sf$in_land <- st_intersects(grid_sf, land, sparse = FALSE) %>% rowSums() > 0  # TRUE for land, FALSE for ocean

# Step 5: Create a land-sea mask based on the land check and ensure correct dimensions
land_sea_mask <- matrix(grid_sf$in_land, nrow = 192, ncol = 96, byrow = FALSE)  # Corrected: 192 longitudes, 96 latitudes
land_sea_mask <- t(land_sea_mask)


### read in IPCC Region information:
#Load reference regions and coastlines:

load("/storage/homefs/ph23v078/Reference_regions/IPCC-WGI-reference-regions-v4_R.rda", verbose = TRUE)

#simplify this object by converting it to a SpatialPolygons class object (i.e., only the polygons are retained and their attributes discarded):
refregions <- as(IPCC_WGI_reference_regions_v4, "SpatialPolygons")

# List of regions to loop over --- excludes ocean basins
regions <- c("GIC", "NWN", "NEN", "WNA", "CNA", "ENA", "NCA", "SCA", "CAR", "NWS",
             "NSA", "NES", "SAM", "SWS", "SES", "SSA", "NEU", "WCE", "EEU", "MED",
             "SAH", "WAF", "CAF", "NEAF", "SEAF", "WSAF", "ESAF", "MDG", "RAR", "WSB",
             "ESB", "RFE", "WCA", "ECA", "TIB", "EAS", "ARP", "SAS", "SEA", "NAU",
             "CAU", "EAU", "SAU", "NZ")
# regions <- c("GIC")



#as before but loop over all regions also, saving everything in an array
path <- "/storage/research/giub_geco/data_2/scratch/phelpap/ModESim"
folders <- list.files(path, full.names = TRUE, pattern = "m[0-9]{3}_tidy") #m001 - m020 in first set

# Extract unique ensemble member identifiers from folder names
ensemble_members <- unique(basename(folders))
# ensemble_members <- "m001_tidy"

### Function to apply land-sea mask
apply_land_mask <- function(grid_data, slm) {
  # Extract the data array from the grid object
  data_array <- grid_data$Data

  # Create a copy of the original data to store masked values
  masked_array <- data_array

  # Loop over the time dimension
  for (t in 1:dim(data_array)[1]) {
    masked_array[t, , ][!land_sea_mask] <- NA
  }

  # Return the modified grid data with the masked values
  grid_data$Data <- masked_array
  return(grid_data)
}

########## for 1420 to 1849, set 1 ###########################

# Initialize a list to store results for all regions
regional_results_1420 <- list()

# Loop over all regions
for (region in regions) {
  # Extract the spatial object corresponding to the current region
  region_object <- refregions[c(region)]

  # Check if the subset is valid
  if (is.null(region_object)) {
    warning(paste("Region not found:", region))
    next
  }

  # Initialize a list to store spatial averages for each file in the current region
  spatial_avg_list <- list()

  for (em in ensemble_members) {
    # Construct the folder path for the current ensemble member
    folder <- file.path(path, em)
    # Construct file paths for the two time periods
    file_1420 <- file.path(folder, "04_result_1420/PCWD_ANNMAX.nc")

    # Check if both files exist
    if (file.exists(file_1420)) {
      # Load grid data for each time period
      grid_1420 <- loadGridData(dataset = file_1420, var = "pcwd_annmax")

      # Apply land-sea mask **before** spatial overlay
      grid_1420 <- apply_land_mask(grid_1420, land_sea_mask)

      # Set spatial projection
      grid_1420 <- setGridProj(grid = grid_1420, proj = proj4string(refregions))

      # Perform spatial overlay
      grid_region_1420 <- overGrid(grid_1420, region_object)

      # Extract data arrays
      data_array_1420 <- grid_region_1420$Data[1:430,,]

      # Compute spatial average for each time step
      spatial_avg <- apply(data_array_1420, 1, function(slice) {
        mean(slice, na.rm = TRUE) # Compute mean for the spatial dimensions, ignoring NA
      })

      # Store the spatial average in the list
      spatial_avg_list[[em]] <- spatial_avg
    } else {
      warning(paste("Missing data for ensemble member:", em))
    }
  }
  # Combine all spatial averages into a single array for the current region
  result_array_region <- do.call(cbind, spatial_avg_list)

  # Store the result for the current region
  regional_results_1420[[region]] <- result_array_region
}
#save calculated list
saveRDS(regional_results_1420, file="~/cwd_global/data/regionalResults_1420_1.RData") #1420_1 for set 1


############## for 1850-2009, set 1 ##########
# Initialize a list to store results for all regions
regional_results_1850 <- list()

# Loop over all regions
for (region in regions) {
  # Extract the spatial object corresponding to the current region
  region_object <- refregions[c(region)]

  # Check if the subset is valid
  if (is.null(region_object)) {
    warning(paste("Region not found:", region))
    next
  }

  # Initialize a list to store spatial averages for each file in the current region
  spatial_avg_list <- list()

  for (em in ensemble_members) {
    # Construct the folder path for the current ensemble member
    folder <- file.path(path, em)
    # Construct file paths for the two time periods
    file_1850 <- file.path(folder, "04_result_1850/PCWD_ANNMAX.nc")

    # Check if both files exist
    if (file.exists(file_1850)) {
      # Load grid data for each time period
      grid_1850 <- loadGridData(dataset = file_1850, var = "pcwd_annmax")

      # Apply land-sea mask **before** spatial overlay
      grid_1850 <- apply_land_mask(grid_1850, land_sea_mask)

      # Set spatial projection
      grid_1850 <- setGridProj(grid = grid_1850, proj = proj4string(refregions))

      # Perform spatial overlay
      grid_region_1850 <- overGrid(grid_1850, region_object)

      # Extract data arrays
      data_array_1850 <- grid_region_1850$Data

      # Compute spatial average for each time step
      spatial_avg <- apply(data_array_1850, 1, function(slice) {
        mean(slice, na.rm = TRUE) # Compute mean for the spatial dimensions, ignoring NA
      })

      # Store the spatial average in the list
      spatial_avg_list[[em]] <- spatial_avg
    } else {
      warning(paste("Missing data for ensemble member:", em))
    }
  }
  # Combine all spatial averages into a single array for the current region
  result_array_region <- do.call(cbind, spatial_avg_list)

  # Store the result for the current region
  regional_results_1850[[region]] <- result_array_region
}

#save calculated list
saveRDS(regional_results_1850, file="~/cwd_global/data/regionalResults_1850_1.RData") #1850_1 for set 1


# # ##################################### Code verification
#
# # Load relevant libraries
# library(dplyr)
# library(rJava)
# library(loadeR.java)
# library(transformeR)
# library(loadeR)
# library(visualizeR)
# library(geoprocessoR)
# library(terra)
# library(ncdf4)
# library(rnaturalearth)
# library(sf)
# library(raster)
# library(ggplot2)
# library(viridis) # For better color scale
# #
# # total_cells <- length(land_sea_mask)
# # ocean_cells <- sum(!land_sea_mask)
# # land_cells <- sum(land_sea_mask)
# #
# # cat("Total Grid Cells:", total_cells, "\n")
# # cat("Land Cells:", land_cells, "(", round(land_cells / total_cells * 100, 2), "% )\n")
# # cat("Ocean Cells:", ocean_cells, "(", round(ocean_cells / total_cells * 100, 2), "% )\n")
# #
#
# # Select one ensemble member for testing
# test_em <- ensemble_members[1]
#
# test_folder <-file.path(path,test_em)
# # Construct file path
# test_file <- file.path(test_folder, "04_result_1420/PCWD_ANNMAX.nc")
#
# ## apply land sea mask on one time step for testing
# # Load data using loadGridData (keeping consistency with your approach)
# grid_test <- loadGridData(dataset = test_file, var = "pcwd_annmax")
# grid_array_test <- grid_test$Data
# test_array <- grid_array_test[1,,]
# test_array[!land_sea_mask] <- NA
#
#
# # Create a copy of the original data to store masked values
# masked_grid_array <- grid_array_test
#
# # Loop over the time dimension
# for (t in 1:dim(grid_array_test)[1]) {
#   masked_grid_array[t, , ][!land_sea_mask] <- NA
# }
#
# # The masked data retains the same 3D structure
# masked_grid_array
#
#
# test_array_2 <- apply_land_mask(grid_test, land_sea_mask)
# test_array_plot <- test_array_2$Data
# plot_data <- test_array_plot[430,,]
#
#
#
#
# r <- raster((plot_data), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
# crs(r) <- "+proj=longlat +datum=WGS84" # Set CRS for raster
# r <- flip(r, direction='y') #Flip the raster to correct orientation
# #Retrieve land polygons and set the same CRS as the raster
# land <- ne_countries(scale = "medium", returnclass = "sf")
# land <- st_transform(land, crs=st_crs(r)) #Ensure CRS alignment
# #Mask the raster with land polygons to remove ocean values
# #r_masked <- mask(r, as(land, "Spatial"))
# #Convert the masked raster to a dataframe for ggplot
# r_df <- as.data.frame(r, xy = TRUE)
# colnames(r_df) <- c("x", "y", "layer")
#
# # Create ggplot
# ggplot() +
#   geom_raster(data = r_df, aes(x = x, y = y, fill = layer), na.rm = TRUE) +
#   scale_fill_viridis_c(option = "C", na.value = "white") +  # Use Viridis color scale
#   labs(title = "GINI Index 1420 epoch (EM)", x = "Longitude", y = "Latitude") +
#   theme_classic() +
#   geom_sf(data = land, fill = NA, color = "black", lwd = 0.5) +  # Overlay coastlines
#   coord_sf(xlim = range(lon), ylim = range(lat), expand = FALSE) +
#   theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
#   guides(fill = guide_colorbar(title = "Gini Index", reverse = FALSE))
