# prerequisites -------------------------------------------------------------------------------

if (FALSE) {
  webshot::install_phantomjs() ## first time only
}

# packages ------------------------------------------------------------------------------------

## data management
library(archive)
library(googledrive)

## data wrangling
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

## GIS packages
library(sf)
library(terra)

## plotting
library(ggplot2)
library(ggspatial)
library(ggtext)
library(scales) ## used for log-scale plotting

# setup ---------------------------------------------------------------------------------------

sf::sf_proj_network(TRUE)
sf::sf_use_s2(TRUE)

extract_mdb <- FALSE ## use TRUE to re-extract from raw data sources (Windows only!)
plot_all <- FALSE ## use TRUE to generate all plots, including exploratory/diagnostic plots
rerun_all <- FALSE ## re-run all analyses, overwriting existing intermediate and output files

## paths
dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
figPath <- "figures" |> fs::dir_create()
outputPath <- "outputs" |> fs::dir_create()

# data download -------------------------------------------------------------------------------

source("01-download-data.R")

# MPB r-value data ----------------------------------------------------------------------------

## mdb extraction only needs to be done once,
## and can only be run on a Windows machine
if (extract_mdb) {
  source("01a-extract-mdb.R")
}

## join all the extracted data into a single table
source("01b-import-mdb-csv.R")

# load pine maps ------------------------------------------------------------------------------

## EOSD (Yemshanov et al. 2012)
yemshanov2012 <- terra::rast(tif_yemshanov2012)

if (plot_all) {
  plot(yemshanov2012)
}

## kNN (Beaudoin et al. 2014)
beaudoin2014 <- terra::rast(tif_beaudoin2014)

if (plot_all) {
  plot(beaudoin2014)
}

## Bleiker 2019
bleiker2019 <- terra::rast(gdb_bleiker2019)

if (plot_all) {
  plot(bleiker2019)
}

## CASFRI

## TODO: use forestData::casfriSpeciesCover() ?

## Coops et al. 2024

## TODO:

## NTEMS (Matasci et al. 2018)

## TODO: use forestData::ntems() ?

## plot them

## TODO: ggplot/cowplot using tidyterra

# load pine introgression (Q) maps  -----------------------------------------------------------

pine_q <- terra::rast(pine_gdb)

if (plot_all) {
  terra::plot(pine_q)
}

# r-values analyses ---------------------------------------------------------------------------

## re-calculate r-values; add BioSIM-generated components
source("02a-Alberta-data-prep.R")

## initial data and model exploration
source("02b-Alberta-explore.R")

##
source("02c-Alberta-analyses.R")

##
source("03-Jasper.R")
