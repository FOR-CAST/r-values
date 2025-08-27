# packages ------------------------------------------------------------------------------------

# library(archive)
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
# geospatial objects for plotting -------------------------------------------------------------

ab_sf <- geodata::gadm("CAN", level = 1, path = dataPath) |>
  sf::st_as_sf() |>
  filter(NAME_1 == "Alberta") |>
  sf::st_geometry()

# pine layers ---------------------------------------------------------------------------------

## TODO: put all data download steps in this script
source("01-data-prep.R")

# MPB r-value data ----------------------------------------------------------------------------

source("01a-extract-mdb.R")

source("01b-import-mdb-csv.R")

# MPB SSI, pine introgression (Q), and MPB winter mortality -----------------------------------

source("02 - Alberta.R")
