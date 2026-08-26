# This script loads IPCC reference regions and adds self-defined additional,
# potentially overlapping regions for analysis
#
# Output: refregion_objects

define_reference_regions <- function() {
  # A) Load IPCC reference regions ----
  load(here::here("data", "Reference_regions", "IPCC-WGI-reference-regions-v4_R.rda"))
  refregions <- as(IPCC_WGI_reference_regions_v4, "SpatialPolygons")

  # List of regions to loop over --- excludes ocean basins
  regions <- c("GIC", "NWN", "NEN", "WNA", "CNA", "ENA", "NCA", "SCA", "CAR", "NWS",
               "NSA", "NES", "SAM", "SWS", "SES", "SSA", "NEU", "WCE", "EEU", "MED",
               "SAH", "WAF", "CAF", "NEAF", "SEAF", "WSAF", "ESAF", "MDG", "RAR", "WSB",
               "ESB", "RFE", "WCA", "ECA", "TIB", "EAS", "ARP", "SAS", "SEA", "NAU",
               "CAU", "EAU", "SAU", "NZ")


  # B) Define additional analysis regions ----

  # define extra_regions on top of IPCC regions for fine grained analysis
  # note that the define regions are allowed to overlap

  # load IPCC regions
  load(here::here("data", "Reference_regions", "IPCC-WGI-reference-regions-v4_R.rda"))


  ## WCE_W and WCE_E ----
  ipcc_sf <- st_as_sf(IPCC_WGI_reference_regions_v4)
  wce_sf  <- ipcc_sf |>
    filter(Acronym == "WCE") |>
    st_make_valid()

  # Example boundary—choose and document the desired longitude.
  wce_split_lon <- 15

  bb <- st_bbox(wce_sf)

  east_box <- st_as_sfc(st_bbox(c(
    xmin = wce_split_lon,
    ymin = unname(bb["ymin"]),
    xmax = unname(bb["xmax"]),
    ymax = unname(bb["ymax"])
  ), crs = st_crs(wce_sf)))

  west_box <- st_as_sfc(st_bbox(c(
    xmin = unname(bb["xmin"]),
    ymin = unname(bb["ymin"]),
    xmax = wce_split_lon,
    ymax = unname(bb["ymax"])
  ), crs = st_crs(wce_sf)))

  wce_e <- st_sf(
    region = "WCE_E",
    geometry = st_intersection(st_geometry(wce_sf), east_box)
  ) |>
    st_make_valid()

  wce_w <- st_sf(
    region = "WCE_W",
    geometry = st_intersection(st_geometry(wce_sf), west_box)
  ) |>
    st_make_valid()

  # drop the added LINESTRING from line boundary
  wce_e <- wce_e[
    st_dimension(wce_e) == 2,
    ,
    drop = FALSE
  ] |> filter(!is.na(region))
  wce_w <- wce_w[
    st_dimension(wce_w) == 2,
    ,
    drop = FALSE
  ] |> filter(!is.na(region))


  ## Switzerland CHE ----
  # che <- ne_countries(
  #   scale = "medium",
  #   country = "Switzerland",
  #   returnclass = "sf"
  # ) |>
  #   st_transform(st_crs(wce_sf)) |>
  #   transmute(region = "CHE")

  ## All European countries ----
  all_european_countries <- ne_countries(
    scale = "medium",
    continent = "Europe",
    returnclass = "sf"
  ) |>
    st_transform(st_crs(wce_sf)) |>
    transmute(region = adm0_iso) |>
    filter(region != "RUS") # remove large country for computational reason

  custom_regions_sf <- bind_rows(wce_e, wce_w, all_european_countries)

  # double check:
  # plot(custom_regions_sf)

  # C) Combined IPCC and self-defined regions ----

  # bind together as a list of SpatialPolygons
  ## transform refregions to list
  refregion_objects <- setNames(
    object =  lapply(regions, function(id) refregions[c(id)]),
    nm = regions
  )
  ## append custom_regions to list
  for (id in custom_regions_sf$region) {
    custom_sf <- custom_regions_sf |> filter(region == id)
    refregion_objects[[id]] <- as(as(custom_sf, "Spatial"), "SpatialPolygons")
  }


  return(refregion_objects)
}

