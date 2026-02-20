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

run_for <- "AB" ## use "NP" to run analyses for national parks; use "AB" to run province-wide
extract_mdb <- FALSE ## use TRUE to re-extract from raw data sources (Windows only!)
plot_all <- FALSE ## use TRUE to generate all plots, including exploratory/diagnostic plots
rerun_all <- FALSE ## re-run all analyses, overwriting existing intermediate and output files

stopifnot(
  run_for %in% c("AB", "NP"),
  is.logical(extract_mdb),
  is.logical(plot_all),
  is.logical(rerun_all)
)

## paths
dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
figPath <- "figures" |> fs::dir_create()
outputPath <- "outputs" |> fs::dir_create()

# data download -------------------------------------------------------------------------------

source("01-download-data.R") ## will prompt for Google authentication

# MPB r-value data ----------------------------------------------------------------------------

## mdb extraction only needs to be done once, and can only be run on a Windows machine

if (run_for == "AB") {
  if (extract_mdb) {
    source("01a-extract-mdb.R")
  }

  ## join all the extracted data into a single table
  source("01b-import-mdb-csv.R")
}

# load pine maps ------------------------------------------------------------------------------

if (run_for == "AB") {
  bleiker2019 <- terra::rast(gdb_bleiker2019)

  if (plot_all) {
    plot(bleiker2019)
  }
}

# load pine introgression (Q) maps  -----------------------------------------------------------

if (run_for == "AB") {
  pine_q <- terra::rast(pine_gdb)

  if (plot_all) {
    terra::plot(pine_q)
  }
}

# r-values analyses ---------------------------------------------------------------------------

if (run_for == "AB") {
  ## re-calculate r-values; add BioSIM-generated components
  source("02a-Alberta-data-prep.R")

  ## initial data and model exploration
  source("02b-Alberta-explore.R")

  ##
  source("02c-Alberta-analyses.R")
} else if (run_for == "NP") {
  ##
  source("03-Jasper.R")
}
