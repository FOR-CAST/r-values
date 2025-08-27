# packages ------------------------------------------------------------------------------------

# library(archive)
# library(cli)
library(dplyr)
# library(elevatr)
# library(fs)
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

# data download -------------------------------------------------------------------------------

source("01-download-data.R")

# MPB r-value data ----------------------------------------------------------------------------

source("01a-extract-mdb.R")

source("01b-import-mdb-csv.R")


# pine layers ---------------------------------------------------------------------------------

source("01c-pine-layers.R")

# MPB SSI, pine introgression (Q), and MPB winter mortality -----------------------------------

source("02-Alberta.R")
