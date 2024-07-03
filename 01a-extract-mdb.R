# packages ------------------------------------------------------------------------------------

library(DBI)
library(odbc)

# setup ---------------------------------------------------------------------------------------

## paths
# cachePath <- "cache"
dataPath <- normalizePath("./data", mustWork = FALSE)
figPath <- "figures"
outputPath <- "outputs"

# if (!dir.exists(cachePath)) dir.create(cachePath)
if (!dir.exists(dataPath)) dir.create(dataPath)
if (!dir.exists(figPath)) dir.create(figPath)
if (!dir.exists(outputPath)) dir.create(outputPath)

# options(reproducible.cachePath = cachePath)

# read in data from MS Access databases -------------------------------------------------------

## connecting to Access databases only works on Windows
stopifnot(.Platform$OS.type == "windows")

mdb_files <- file.path(dataPath, "AB", "mdb") |>
  list.files(full.names = TRUE, pattern = "[.]mdb$", recursive = TRUE) |>
  grep("/Copy of", x = _, invert = TRUE, value = TRUE) |> ## filter unwanted files
  grep("rollup", x = _, invert = TRUE, value = TRUE) |>   ## filter unwanted files
  normalizePath()

purrr::walk(mdb_files, function(mdb) {
  tmp_mdb <- tempfile(fileext = ".mdb")
  file.copy(mdb, tmp_mdb) ## work on a copy; needed with network drive
  connection_string <- paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=", tmp_mdb)
  con <- dbConnect(odbc::odbc(), .connection_string = connection_string, timeout = 10)

  ## read in the table as an R data.frame and export to .csv
  db_tbls <- dbListTables(con)

  if ("mpb_trees" %in% db_tbls) {
    mpb_trees <- dbReadTable(con, "mpb_trees")
    write.csv(mpb_trees, paste0(tools::file_path_sans_ext(mdb), "_mpb_trees.csv"), row.names = FALSE)
  }

  if ("mpb_survey_info" %in% db_tbls) {
    surv_info <- dbReadTable(con, "mpb_survey_info")
    write.csv(surv_info, paste0(tools::file_path_sans_ext(mdb), "_mpb_survey_info.csv"), row.names = FALSE)
  }

  if ("mpb_site" %in% db_tbls) {
    mpb_site <- dbReadTable(con, "mpb_site")
    write.csv(mpb_site, paste0(tools::file_path_sans_ext(mdb), "_mpb_site.csv"), row.names = FALSE)
  }

  ## cleanup
  dbDisconnect(con)
  unlink(tmp_mdb)
})

csv_files <- c(
  paste0(tools::file_path_sans_ext(mdb_files), "_mpb_trees.csv"),
  paste0(tools::file_path_sans_ext(mdb_files), "_mpb_survey_info.csv"),
  paste0(tools::file_path_sans_ext(mdb_files), "_mpb_site.csv")
)
csv_files <- csv_files[file.exists(csv_files)] ## not all tables/files exist, so filter them out
csv_files <- fs::path_rel(csv_files, dataPath)

withr::with_dir(dataPath, archive::archive_write_files("extracted_mdb_tables.zip", csv_files))

googledrive::drive_put(file.path(dataPath, "extracted_mdb_tables.zip"), drive_id)
## uplooaded to '16OK48u6g5-JdWvEFhIJ0vMxvcK-k6z5l'
