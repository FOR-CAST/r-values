# packages ------------------------------------------------------------------------------------

library(archive)
library(geodata)
library(googledrive)
library(purrr)
library(sf)
library(terra)

# setup ---------------------------------------------------------------------------------------

## paths
dataPath <- normalizePath("./data", mustWork = FALSE)
figPath <- "figures"
outputPath <- "outputs"

if (!dir.exists(dataPath)) dir.create(dataPath)
if (!dir.exists(figPath)) dir.create(figPath)
if (!dir.exists(outputPath)) dir.create(outputPath)

## set map projection
latlon <- crs("epsg:4326")
targetCRS <- crs(paste(
  "+proj=aea +lat_1=49 +lat_2=67 +lat_0=0 +lon_0=-112 +x_0=0 +y_0=0",
  "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
))

# get data from google drive ------------------------------------------------------------------

drive_id <- as_id("1EiproEknMuuze5c6U_1tCptB8UM4z5YI")

## create local data directory structure to match that on Google Drive
drive_dirs <- drive_ls(drive_id, recursive = FALSE, type = "folder")

purrr::walk(drive_dirs$name, function(d) {
  if (!dir.exists(file.path(dataPath, d))) {
    dir.create(file.path(dataPath, d), recursive = TRUE)
  }
})

## fetch file info for ALL files; will refer to this data frame as needed to download files
drive_files <- drive_ls(drive_id, recursive = TRUE)

## TODO: download files

# Alberta administrative boundaries -----------------------------------------------------------

## Canadian provincial/territorial boundaries
can1.latlon <- geodata::gadm(country = "CAN", level = 1, path = dataPath, version = "4.1") |>
  st_as_sf()
can1 <- st_transform(can1.latlon, targetCRS)

ab.latlon <- can1.latlon[can1.latlon$NAME_1 == "Alberta", ]
ab <- can1[can1$NAME_1 == "Alberta", ]

# get pine maps -------------------------------------------------------------------------------

## TODO: issue #3

# get SLR values from BIOSIM ------------------------------------------------------------------

## TODO: issue #2

