# packages ------------------------------------------------------------------------------------

library(archive)
library(geodata)
library(googledrive)
library(purrr)
library(sf)
library(terra)

library(reproducible)
library(LandR)

## NOTE: several suggested packages used here but renv complains:
#
# package  installed recorded used
# BioSIM   y         y        n
# cluster  y         y        n
# foreign  y         y        n
# httpuv   y         y        n
# J4R      y         y        n
# later    y         y        n
# promises y         y        n
# rpart    y         y        n
# spatial  y         y        n
# XML      y         y        n

# setup ---------------------------------------------------------------------------------------

## paths
cachePath <- "cache"
dataPath <- normalizePath("./data", mustWork = FALSE)
figPath <- "figures"
outputPath <- "outputs"

if (!dir.exists(cachePath)) dir.create(cachePath)
if (!dir.exists(dataPath)) dir.create(dataPath)
if (!dir.exists(figPath)) dir.create(figPath)
if (!dir.exists(outputPath)) dir.create(outputPath)

options(reproducible.cachePath = cachePath)

## set map projection
latlon <- crs("epsg:4326")
targetCRS <- crs(paste(
  "+proj=aea +lat_1=49 +lat_2=67 +lat_0=0 +lon_0=-112 +x_0=0 +y_0=0",
  "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
))

source("R/helpers.R")

# get data from google drive ------------------------------------------------------------------

drive_id <- as_id("1EiproEknMuuze5c6U_1tCptB8UM4z5YI")

## NOTE: use overwrite=TRUE if e.g., data updated on Google Drive
all_drive_files <- workflowtools::drive_download_folder(drive_id, dataPath, overwrite = FALSE)

# Alberta administrative boundaries -----------------------------------------------------------

## Canadian provincial/territorial boundaries
can1.latlon <- geodata::gadm(country = "CAN", level = 1, path = dataPath, version = "4.1") |>
  st_as_sf()
can1 <- st_transform(can1.latlon, targetCRS)

ab.latlon <- can1.latlon[can1.latlon$NAME_1 == "Alberta", ]
ab <- can1[can1$NAME_1 == "Alberta", ]

# create DEM for use with BioSIM --------------------------------------------------------------

## follows approach taken in LandR_MPB_studyArea and mpbClimateData modules
absk <- subset(can1.latlon, NAME_1 %in% c("Alberta", "Saskatchewan")) |>
  sf::st_transform(targetCRS)

studyAreaReporting <- mpbutils::mpbStudyArea(
  ecoregions = c(112, 120, 122, 124, 126),
  targetCRS = targetCRS,
  cPath = cachePath,
  dPath = dataPath
) |>
  sf::st_intersection(absk) |>
  sf::st_union()

rasterToMatchReporting <- Cache(
  LandR::prepInputsLCC,
  year = 2005, ## TODO: use 2010?
  studyArea = studyAreaReporting,
  destinationPath = dataPath,
  filename2 = NULL
)

DEM <- Cache(
  LandR::prepInputsCanDEM,
  rasterToMatch = rasterToMatchReporting,
  studyArea = studyAreaReporting,
  destinationPath = dataPath
) |>
  writeRaster(file.path(outputPath, "DEM_ABSK_studyArea.tif")) ## TODO: different file format?

# get pine maps -------------------------------------------------------------------------------

## TODO: issue #3

## EOSD (Yemshanov et al. 2012)
url_yemshanov2012 <- as_id("11g02bDnEt6U_xXtcLWLmqzWLleR_c54F")
yemshanov2012 <- prepInputs(
  url = url_yemshanov2012,
  destinationPath = dataPath,
  fun = "terra::rast",
  cropTo = ab,
  maskTo = ab,
  targetFile = "Yemshanov_pine_map.flt"
) |>
  Cache() |>
  writeRaster(file.path(dataPath, "Yemshanov_pine_map.tif"), overwrite = TRUE)
crs_yemshanov2012 <- crs(yemshanov2012)

## kNN (Beaudoin et al. 2014)
sppEquiv <- LandR::sppEquivalencies_CA
sppEquiv <- sppEquiv[KNN %in% c("Pinu_Ban", "Pinu_Con"), ]

beaudoin2014 <- prepSpeciesLayers_KNN(
  destinationPath = dataPath,
  outputPath = dataPath,
  url = NULL,
  studyArea = ab,
  rasterToMatch = yemshanov2012,
  sppEquiv = sppEquiv,
  sppEquivCol = "KNN",
  thresh = 10, ## i.e., minimum 10% cover
  year = 2011
) |>
  Cache() |>
  writeRaster(file.path(dataPath, "Beaudoin_pine_map.tif"), overwrite = TRUE)
crs_beaudoin2014 <- crs(beaudoin2014)

## Bleiker 2019
url_bleiker2019 <- as_id("15EzncjIR_dn5v6hruoVbsUQVF706nTEL")
bleiker2019 <- prepInputs_ABPine(
  url = url_bleiker2019,
  destinationPath = dataPath,
  layerNames = "OVERSTOREY_PINE",
  rasterToMatch = yemshanov2012 ## TODO: could do better than 250m
) |>
  Cache() |>
  writeRaster(file.path(dataPath, "Bleiker_pine_map.tif"), overwrite = TRUE)
crs_bleiker2019 <- crs(bleiker2019)

## CASFRI

## TODO: use forestData::casfriSpeciesCover() ?

## Coops et al. 2024

## TODO:

## NTEMS (Matasci et al. 2018)

## TODO: use forestData::ntems() ?

## plot them

## TODO: ggplot/cowplot using tidyterra

# get SLR values from BIOSIM ------------------------------------------------------------------

## TODO: issue #2
