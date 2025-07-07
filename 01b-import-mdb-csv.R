# packages ------------------------------------------------------------------------------------

# library(archive)
library(dplyr)
# library(fs)
# library(googledrive)
# library(purrr)

# setup ---------------------------------------------------------------------------------------

## paths
dataPath <- normalizePath("./data", mustWork = FALSE)
figPath <- "figures"
outputPath <- "outputs"

dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
figPath <- "figures" |> fs::dir_create()
outputPath <- "outputs" |> fs::dir_create()

# validate and merge csv tables ---------------------------------------------------------------

csv_zip <- file.path(dataPath, "extracted_mdb_tables.zip")

if (!file.exists(csv_zip)) {
  googledrive::drive_download(googledrive::as_id("16OK48u6g5-JdWvEFhIJ0vMxvcK-k6z5l"))
  archive::archive_extract(csv_zip, dataPath)
}

csv_files <- dataPath |>
  list.files(pattern = "_mpb_(site|trees)[.]csv", full.names = TRUE, recursive = TRUE) |>
  fs::path_rel()

## read in each (set of) files, validate, and merge
dirname(csv_files) |>
  unique() |>
  purrr::walk(function(d) {
    ## if mpb_site table missing, use mpb_survey_info
    survey_site_csv <- file.path(d, list.files(d, pattern = "_mpb_site[.]csv"))
    if (length(survey_site_csv) == 0) {
      survey_site_csv <- file.path(d, list.files(d, pattern = "_mpb_survey_info[.]csv"))
    }

    mpb_trees_csv <- file.path(d, list.files(d, pattern = "_mpb_trees[.]csv"))

    purrr::walk2(survey_site_csv, mpb_trees_csv, function(fsite, ftrees) {
      ## check the csvs are correctly paired before joining
      stopifnot(identical(sub("(mpb_site|mpb_survey_info)", "mpb_trees", basename(fsite)), basename(ftrees)))

      survey_site <- read.csv(fsite) |>
        mutate(infestation = as.character(infestation))

      mpb_trees <- read.csv(ftrees) |>
        mutate(infestation = as.character(infestation))

      mpb_site_trees <- left_join(mpb_trees, survey_site, by = "siteID") |>
        rename(infestation = infestation.y, site_nbr = site_nbr.y) |>
        mutate(infestation.x = NULL, site_nbr.x = NULL)

      ## TODO: more validation
      ## sample_tree_nbr and tree_nbr

      sub("(mpb_site|mpb_survey_info)", "mpb_site_trees_cleaned", fsite) |>
        write.csv(mpb_site_trees, file = _, row.names = FALSE)
    })
    ## TODO: merge tables
  })
