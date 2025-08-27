# get data from google drive ------------------------------------------------------------------

drive_id <- googledrive::as_id("1EiproEknMuuze5c6U_1tCptB8UM4z5YI")

if (FALSE) {
  ## look up file ids etc. as needed
  all_drive_files <- googledrive::drive_ls(drive_id, recursive = TRUE)
}

## NOTE: use overwrite=TRUE if e.g., data updated on Google Drive
downloaded_files <- workflowtools::drive_download_folder(
  drive_id,
  dataPath,
  overwrite = FALSE
)

## e.g., to re-download a subdirectory only:
# withr::with_dir(
#   file.path(dataPath, "Brett"),
#   workflowtools::drive_download_folder(
#     as_id("12aHAFqjL40ly6x9tvfjj_96HE0ZjOQyl"), ## get drive folder ID from the 'share' link/url
#     file.path(dataPath, "Brett"),
#     overwrite = TRUE
#   )
# )

## ensure the 2011 beetle year data are present locally
## (they were uploaded after the initial drive upload)
f_zip_2011 <- file.path(dataPath, "AB", "2011_Population_forecast_r-value.zip")
d_zip_2011 <- file.path(dataPath, "AB", "mdb", "SourceData2009to2011") ## dest dir
if (!file.exists(f_zip_2011)) {
  googledrive::as_id("1z2KmKFAar-G0-5iJGdi1uhcv5gSftfPu") |>
    googledrive::drive_download(f_zip_2011)
}

if (!dir.exists(file.path(d_zip_2011, "Population forecast (r value)"))) {
  archive::archive_extract(archive = f_zip_2011, dir = d_zip_2011)
}

# geospatial objects for plotting -------------------------------------------------------------

ab_sf <- try(
  geodata::gadm("CAN", level = 1, path = dataPath) |>
    sf::st_as_sf() |>
    filter(NAME_1 == "Alberta") |>
    sf::st_geometry()
)

if (inherits(ab_sf, "try-error")) {
  gadm_can_rds <- file.path(dataPath, "gadm", "gadm41_CAN_1_pk.rds")
  if (!file.exists(gadm_can_rds)) {
    as_id("1zTMd5p9jufwRVGkD2IBjeLu20nFj6MsS") |>
      drive_download(path = gadm_can_rds)
  }

  ab_sf <- readRDS(gadm_can_rds) |>
    sf::st_as_sf() |>
    filter(NAME_1 == "Alberta") |>
    sf::st_geometry()
}
