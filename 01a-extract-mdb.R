# packages ------------------------------------------------------------------------------------

library(DBI)
library(googledrive)
library(odbc)

# setup ---------------------------------------------------------------------------------------

## paths
dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
figPath <- "figures" |> fs::dir_create()
outputPath <- "outputs" |> fs::dir_create()

# read in data from MS Access databases -------------------------------------------------------

## connecting to Access databases only works on Windows
stopifnot(.Platform$OS.type == "windows")

mdb_files <- file.path(dataPath, "AB", "mdb") |>
  list.files(full.names = TRUE, pattern = "[.](accdb|mdb)$", recursive = TRUE) |>
  ## filter unwanted files:
  grep("/Copy of", x = _, invert = TRUE, value = TRUE) |>
  grep("rollup", x = _, invert = TRUE, value = TRUE) |>
  grep("/source/", x = _, invert = TRUE, value = TRUE) |>
  grep("/PopForecast_Master2009[.]mdb", x = _, invert = TRUE, value = TRUE) |>
  grep("mpb_survey_Slave_2009_woSP", x = _, invert = TRUE, value = TRUE) |>
  normalizePath()

## cleanup any pre-existing files from previous runs
dirname(mdb_files) |>
  unique() |>
  vapply(function(d) fs::path(d, c("source", "site", "tree")), character(3)) |>
  purrr::walk(function(d) if (dir.exists(d)) unlink(d, recursive = TRUE))

## extract mdb tables and output to csv
purrr::walk(mdb_files, function(mdb) {
  tmp_mdb <- tempfile(fileext = paste0(".", tools::file_ext(mdb)))
  file.copy(mdb, tmp_mdb) ## work on a copy; needed with network drive
  connection_string <- paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=", tmp_mdb)
  con <- dbConnect(odbc::odbc(), .connection_string = connection_string, timeout = 10)

  ## read in the table as an R data.frame and export to .csv
  db_tbls <- dbListTables(con)

  mdb_dir <- dirname(mdb)
  out_dir <- fs::dir_create(outputPath, fs::path_rel(mdb_dir, dataPath))

  ## organize the accdb/mdb files into a 'source/' subdirectory;
  ## organize the site-level .csv files into a 'site/' subdirectory;
  ## organize the tree-level .csv files into a 'tree/' subdirectory;
  fs::dir_create(out_dir, "source")
  fs::dir_create(out_dir, "site")
  fs::dir_create(out_dir, "tree")

  file.copy(mdb, file.path(out_dir, "source", basename(mdb)))

  if ("mpb_trees" %in% db_tbls) {
    mpb_trees <- dbReadTable(con, "mpb_trees")
    csv_file <- file.path(out_dir, "tree", paste0(tools::file_path_sans_ext(basename(mdb)), "_mpb_trees.csv")) |>
      normalizePath(mustWork = FALSE)
    write.csv(mpb_trees, csv_file, row.names = FALSE)
  }

  ## site and survey_info are equivalent, so put in the 'site/' subdir
  if ("mpb_survey_info" %in% db_tbls) {
    surv_info <- dbReadTable(con, "mpb_survey_info")
    csv_file <- file.path(out_dir, "site", paste0(tools::file_path_sans_ext(basename(mdb)), "_mpb_survey_info.csv")) |>
      normalizePath(mustWork = FALSE)
    write.csv(surv_info, csv_file, row.names = FALSE)
  }

  if ("mpb_site" %in% db_tbls) {
    mpb_site <- dbReadTable(con, "mpb_site")
    csv_file <- file.path(out_dir, "site", paste0(tools::file_path_sans_ext(basename(mdb)), "_mpb_site.csv")) |>
      normalizePath(mustWork = FALSE)
    write.csv(mpb_site, csv_file, row.names = FALSE)
  }

  ## cleanup
  dbDisconnect(con)
  unlink(tmp_mdb)
})

## expected outputs from previous
csv_files <- c(
  file.path(dirname(mdb_files), "tree", paste0(tools::file_path_sans_ext(basename(mdb_files)), "_mpb_trees.csv")),
  file.path(dirname(mdb_files), "site", paste0(tools::file_path_sans_ext(basename(mdb_files)), "_mpb_survey_info.csv")),
  file.path(dirname(mdb_files), "site", paste0(tools::file_path_sans_ext(basename(mdb_files)), "_mpb_site.csv"))
)
csv_files <- fs::path(outputPath, fs::path_rel(csv_files, dataPath))
csv_files <- csv_files[file.exists(csv_files)] ## not all tables/files exist, so filter them out

## archive and upload ---------------------------------------------------------

## remove any previously created versions
prev_archive <- file.path(dataPath, "extracted_mdb_tables.zip")
if (file.exists(prev_archive)) {
  file.remove(prev_archive)
}

## create archive
withr::with_dir(outputPath, {
  archive::archive_write_files("extracted_mdb_tables.zip", fs::path_rel(csv_files, outputPath))
})

## upload to Google Drive
drive_id <- as_id("16OK48u6g5-JdWvEFhIJ0vMxvcK-k6z5l")
drive_update(file = drive_id, media = file.path(outputPath, "extracted_mdb_tables.zip"))
