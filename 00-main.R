# packages ------------------------------------------------------------------------------------

## data management
library(archive)
library(googledrive)

## data wrangling
library(dplyr)
library(purrr)
library(tidyr)

## GIS packages
library(sf)
library(terra)

## plotting
library(ggplot2)
library(ggspatial)

# setup ---------------------------------------------------------------------------------------

## paths
dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
figPath <- "figures" |> fs::dir_create()
outputPath <- "outputs" |> fs::dir_create()

# data download -------------------------------------------------------------------------------

source("01-download-data.R")

# MPB r-value data ----------------------------------------------------------------------------

## mdb extraction only needs to be done once,
## and can only be run on a Windows machine
source("01a-extract-mdb.R")

## join all the extracted data into a single table
source("01b-import-mdb-csv.R")

# pine layers ---------------------------------------------------------------------------------

source("01c-pine-layers.R")

# MPB SSI, pine introgression (Q), and MPB winter mortality -----------------------------------

source("02a-Alberta-data-prep.R")
source("02b-Alberta-explore.R")
source("02c-Alberta-analyses.R")
