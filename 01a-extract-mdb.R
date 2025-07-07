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
  grep("/Copy of", x = _, invert = TRUE, value = TRUE) |> ## filter unwanted files
  grep("rollup", x = _, invert = TRUE, value = TRUE) |>   ## filter unwanted files
  normalizePath()

purrr::walk(mdb_files, function(mdb) {
  tmp_mdb <- tempfile(fileext = paste0(".", tools::file_ext(mdb)))
  file.copy(mdb, tmp_mdb) ## work on a copy; needed with network drive
  connection_string <- paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=", tmp_mdb)
  con <- dbConnect(odbc::odbc(), .connection_string = connection_string, timeout = 10)

  ## read in the table as an R data.frame and export to .csv
  db_tbls <- dbListTables(con)

  ## organize the accdb/mdb files into a 'source/' subdirectory;
  ## organize the site-level .csv files into a 'site/' subdirectory;
  ## organize the tree-level .csv files into a 'tree/' subdirectory;
  fs::dir_create(dirname(mdb), "source")
  fs::dir_create(dirname(mdb), "site")
  fs::dir_create(dirname(mdb), "tree")

  file.copy(mdb, file.path(dirname(mdb), "source", basename(mdb)))

  if ("mpb_trees" %in% db_tbls) {
    mpb_trees <- dbReadTable(con, "mpb_trees")
    csv_file <- file.path(dirname(mdb), "tree", paste0(tools::file_path_sans_ext(basename(mdb)), "_mpb_trees.csv")) |>
      normalizePath(mustWork = FALSE)
    write.csv(mpb_trees, csv_file, row.names = FALSE)
  }

  ## site and survey_info are equivalent, so put in the 'site/' subdir
  if ("mpb_survey_info" %in% db_tbls) {
    surv_info <- dbReadTable(con, "mpb_survey_info")
    csv_file <- file.path(dirname(mdb), "site", paste0(tools::file_path_sans_ext(basename(mdb)), "_mpb_survey_info.csv")) |>
      normalizePath(mustWork = FALSE)
    write.csv(surv_info, csv_file, row.names = FALSE)
  }

  if ("mpb_site" %in% db_tbls) {
    mpb_site <- dbReadTable(con, "mpb_site")
    csv_file <- file.path(dirname(mdb), "site", paste0(tools::file_path_sans_ext(basename(mdb)), "_mpb_site.csv")) |>
      normalizePath(mustWork = FALSE)
    write.csv(mpb_site, csv_file, row.names = FALSE)
  }

  ## cleanup
  dbDisconnect(con)
  unlink(tmp_mdb)
})

csv_files <- c(
  file.path(dirname(mdb_files), "tree", paste0(tools::file_path_sans_ext(basename(mdb_files)), "_mpb_trees.csv")),
  file.path(dirname(mdb_files), "site", paste0(tools::file_path_sans_ext(basename(mdb_files)), "_mpb_survey_info.csv")),
  file.path(dirname(mdb_files), "site", paste0(tools::file_path_sans_ext(basename(mdb_files)), "_mpb_site.csv"))
)
csv_files <- csv_files[file.exists(csv_files)] ## not all tables/files exist, so filter them out
csv_files <- fs::path_rel(csv_files, dataPath)

withr::with_dir(dataPath, {
  archive::archive_write_files("extracted_mdb_tables.zip", csv_files)
})

drive_id <- as_id("16OK48u6g5-JdWvEFhIJ0vMxvcK-k6z5l")
drive_update(file = drive_id, media = file.path(dataPath, "extracted_mdb_tables.zip"))
