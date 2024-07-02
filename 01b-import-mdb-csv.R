# packages ------------------------------------------------------------------------------------

library(archive)
library(dplyr)
library(googledrive)

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


# validate and merge csv tables ---------------------------------------------------------------

csv_zip <- file.path(dataPath, "extracted_mdb_tables.zip")

if (!file.exists(csv_zip)) {
  drive_download(as_id("16OK48u6g5-JdWvEFhIJ0vMxvcK-k6z5l"))
  archive_extract(csv_zip, dataPath)
}

csv_files <- list.files(dataPath, pattern = "_(mpb_trees|mpb_survey_info|mpb_site)[.]csv", recursive = TRUE)
csv_files <- file.path(dataPath, csv_files)

## TODO: read in each (set of) files, validate, and merge
