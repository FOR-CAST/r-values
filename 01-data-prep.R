# packages ------------------------------------------------------------------------------------

library(archive)
library(dplyr)
# library(elevatr)
# library(geodata)
library(ggplot2)
library(googledrive)
library(purrr)
# library(RCurl)
library(sf)
library(terra)
# library(XML)

# setup ---------------------------------------------------------------------------------------

## paths
dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
figPath <- "figures" |> fs::dir_create()
outputPath <- "outputs" |> fs::dir_create()

## set map projection
latlon <- crs("epsg:4326")
targetCRS <- crs(paste(
  "+proj=aea +lat_1=49 +lat_2=67 +lat_0=0 +lon_0=-112 +x_0=0 +y_0=0",
  "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
))

# get data from google drive ------------------------------------------------------------------

drive_id <- as_id("1EiproEknMuuze5c6U_1tCptB8UM4z5YI")

all_drive_files <- drive_ls(drive_id, recursive = TRUE) ## can be used to look up file ids etc.

## NOTE: use overwrite=TRUE if e.g., data updated on Google Drive
downloaded_files <- workflowtools::drive_download_folder(
  drive_id,
  dataPath,
  overwrite = FALSE
)

## e.g., to re-download a subdirectory only:
# withr::with_dir(
#   file.path(dataPath, "Brett"),
#   workflowtools::drive_download_folder(
#     as_id("12aHAFqjL40ly6x9tvfjj_96HE0ZjOQyl"), ## get drive folder ID from the 'share' link/url
#     file.path(dataPath, "Brett"),
#     overwrite = TRUE
#   )
# )

## ensure the 2011 beetle year data are present locally
## (they were uploaded after the initial drive upload)
f_zip_2011 <- file.path(dataPath, "AB", "2011_Population_forecast_r-value.zip")
d_zip_2011 <- file.path(dataPath, "AB", "mdb", "SourceData2009to2011") ## dest dir
if (!file.exists(f_zip_2011)) {
  as_id("1z2KmKFAar-G0-5iJGdi1uhcv5gSftfPu") |>
    drive_download(f_zip_2011)
}

if (!dir.exists(file.path(d_zip_2011, "Population forecast (r value)"))) {
  archive::archive_extract(
    archive = f_zip_2011,
    dir = d_zip_2011
  )
}

# geospatial objects for plotting -------------------------------------------------------------

ab_sf <- geodata::gadm("CAN", level = 1, path = dataPath) |>
  sf::st_as_sf() |>
  filter(NAME_1 == "Alberta") |>
  sf::st_geometry()

# get pine maps -------------------------------------------------------------------------------

## TODO: issue #3

## EOSD (Yemshanov et al. 2012)
url_yemshanov2012 <- as_id("11g02bDnEt6U_xXtcLWLmqzWLleR_c54F")
yemshanov2012 <- prepInputs(
  url = url_yemshanov2012,
  destinationPath = dataPath,
  fun = "terra::rast",
  cropTo = ab_sf,
  maskTo = ab_sf,
  targetFile = "Yemshanov_pine_map.flt"
) |>
  Cache() |>
  writeRaster(file.path(dataPath, "Yemshanov_pine_map.tif"), overwrite = TRUE)
crs_yemshanov2012 <- crs(yemshanov2012)

## kNN (Beaudoin et al. 2014)
tif_beaudoin2014 <- file.path(dataPath, "Beaudoin_pine_map.tif")

if (file.exists(tif_beaudoin2014)) {
  beaudoin2014 <- terra::rast(tif_beaudoin2014)
} else {
  url_beaudoin2014 <- paste0(
    "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
    "canada-forests-attributes_attributs-forests-canada/2011",
    "-attributes_attributs-2011/"
  )

  fileURLs <- RCurl::getURL(
    url_beaudoin2014,
    dirlistonly = TRUE,
    .opts = list(followlocation = TRUE)
  )
  fileNames <- XML::getHTMLLinks(fileURLs)
  fileNames <- grep("(Species_Pinu_Ban|Species_Pinu_Con)_.*\\.tif", fileNames, value = TRUE)
  utils::download.file(
    url = paste0(url_beaudoin2014, fileNames),
    destfile = file.path(dataPath, fileNames)
  )

  beaudoin2014 <- file.path(dataPath, fileNames) |>
    grep("[.]tif$", x = _, value = TRUE) |>
    terra::rast()
  beaudoin2014 <- terra::crop(
    x = beaudoin2014,
    y = sf::st_transform(ab_sf, terra::crs(beaudoin2014))
  )
  beaudoin2014 <- terra::mask(
    x = beaudoin2014,
    mask = sf::st_transform(ab_sf, terra::crs(beaudoin2014)) |> terra::vect()
  )
  terra::set.names(beaudoin2014, c("Pinu_Ban", "Pinu_Con"))

  terra::plot(beaudoin2014)

  beaudoin2014 <- terra::writeRaster(beaudoin2014, tif_beaudoin2014)
}

## Bleiker 2019
gdb_bleiker2019 <- file.path(dataPath, "AB_PineVolumes_Lambert.gdb")
zip_bleiker2019 <- paste0(gdb_bleiker2019, ".zip")

if (!file.exists(zip_bleiker2019)) {
  as_id("15EzncjIR_dn5v6hruoVbsUQVF706nTEL") |>
    drive_download(path = zip_bleiker2019, overwrite = TRUE)
}

if (!(file.exists(gdb_bleiker2019) || dir.exists(gdb_bleiker2019))) {
  archive::archive_extract(zip_bleiker2019, dataPath)
}

bleiker2019 <- sf::st_read(gdb_bleiker2019, layer = "OVERSTOREY_PINE") |>
  sf::st_make_valid()

## CASFRI

## TODO: use forestData::casfriSpeciesCover() ?

## Coops et al. 2024

## TODO:

## NTEMS (Matasci et al. 2018)

## TODO: use forestData::ntems() ?

## plot them

## TODO: ggplot/cowplot using tidyterra
