# packages ------------------------------------------------------------------------------------

library(DBI)
library(odbc)
library(ggplot2)
library(ggspatial)

# create DEM for use with BioSIM --------------------------------------------------------------

## follows approach taken in LandR_MPB_studyArea and mpbClimateData modules
absk <- geodata::gadm(country = "CAN", level = 1, path = dataPath) |>
  sf::st_as_sf() |>
  subset(x = _, NAME_1 %in% c("Alberta", "Saskatchewan")) |>
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


# read in data from MS Access databases -------------------------------------------------------

## TODO
