# get data from google drive ------------------------------------------------------------------

local({
  drive_id <- googledrive::as_id("1EiproEknMuuze5c6U_1tCptB8UM4z5YI")

  drive_dirs <- googledrive::drive_ls(drive_id, recursive = FALSE, type = "folder")

  ## NOTE: use overwrite=TRUE if e.g., data updated on Google Drive
  apply(drive_dirs, 1, function(d) {
    dPath <- file.path(dataPath, d$name)
    if (!dir.exists(dPath)) {
      workflowtools::drive_download_folder(
        drive_folder = d$id,
        path = dPath,
        overwrite = FALSE
      )
    } else {
      NULL
    }
  })

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
})

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
    googledrive::as_id("1zTMd5p9jufwRVGkD2IBjeLu20nFj6MsS") |>
      googledrive::drive_download(path = gadm_can_rds)
  }

  ab_sf <- readRDS(gadm_can_rds) |>
    sf::st_as_sf() |>
    filter(NAME_1 == "Alberta") |>
    sf::st_geometry()

  rm(gadm_can_rds)
}

# pine maps -----------------------------------------------------------------------------------

## TODO: issue #3

tif_yemshanov2012 <- file.path(dataPath, "Yemshanov_pine_map.tif")
tif_beaudoin2014 <- file.path(dataPath, "Beaudoin_pine_map.tif")
gdb_bleiker2019 <- file.path(dataPath, "AB_PineVolumes_Lambert.gdb")

local({
  ## EOSD (Yemshanov et al. 2012)
  zip_yemshanov2012 <- file.path(dataPath, "Yemshanov_pine_map.zip")
  if (!file.exists(zip_yemshanov2012)) {
    googledrive::as_id("11g02bDnEt6U_xXtcLWLmqzWLleR_c54F") |>
      googledrive::drive_download(path = zip_yemshanov2012, overwrite = TRUE)
  }

  if (!file.exists(tif_yemshanov2012)) {
    archive::archive_extract(zip_yemshanov2012, dataPath)

    yemshanov2012 <- file.path(dataPath, "Yemshanov_pine_map", "Yemshanov_pine_map.flt") |>
      terra::rast()

    yemshanov2012 <- terra::crop(
      x = yemshanov2012,
      y = terra::project(terra::vect(ab_sf), yemshanov2012),
      mask = TRUE
    ) |>
      terra::writeRaster(tif_yemshanov2012, overwrite = TRUE)
  }

  ## kNN (Beaudoin et al. 2014)
  if (!file.exists(tif_beaudoin2014)) {
    url_beaudoin2014 <- paste0(
      "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
      "canada-forests-attributes_attributs-forests-canada/2011",
      "-attributes_attributs-2011/"
    )

    fileURLs <- RCurl::getURL(
      url_beaudoin2014,
      dirlistonly = TRUE,
      .opts = list(followlocation = TRUE)
    )
    fileNames <- XML::getHTMLLinks(fileURLs)
    fileNames <- grep("(Species_Pinu_Ban|Species_Pinu_Con)_.*\\.tif", fileNames, value = TRUE)
    utils::download.file(
      url = paste0(url_beaudoin2014, fileNames),
      destfile = file.path(dataPath, fileNames)
    )

    beaudoin2014 <- file.path(dataPath, fileNames) |>
      grep("[.]tif$", x = _, value = TRUE) |>
      terra::rast()
    beaudoin2014 <- terra::crop(
      x = beaudoin2014,
      y = sf::st_transform(ab_sf, terra::crs(beaudoin2014))
    )
    beaudoin2014 <- terra::mask(
      x = beaudoin2014,
      mask = sf::st_transform(ab_sf, terra::crs(beaudoin2014)) |> terra::vect()
    )
    terra::set.names(beaudoin2014, c("Pinu_Ban", "Pinu_Con"))

    beaudoin2014 <- terra::writeRaster(beaudoin2014, tif_beaudoin2014)
  }

  ## Bleiker 2019
  zip_bleiker2019 <- paste0(gdb_bleiker2019, ".zip")

  if (!file.exists(zip_bleiker2019)) {
    googledrive::as_id("15EzncjIR_dn5v6hruoVbsUQVF706nTEL") |>
      googledrive::drive_download(path = zip_bleiker2019, overwrite = TRUE)
  }

  if (!(file.exists(gdb_bleiker2019) || dir.exists(gdb_bleiker2019))) {
    archive::archive_extract(zip_bleiker2019, dataPath)
  }
})

# pine introgression (Q) maps -----------------------------------------------------------------

q_map <- file.path(dataPath, "HybridPrediction_1000m.lpkx") ## download to here
q_map_dir <- file.path(dataPath, fs::path_ext_remove(basename(q_map))) ## extract to here

if (!file.exists(q_map)) {
  googledrive::as_id("1L2UPbWoXH5_uXuRkaKBU9nZl5ZUvZCnW") |>
    googledrive::drive_download(path = q_map)

  ## .lpkx is a 7zip folder that bundles .gdb files, so we can simply extract:
  ##  - the p20 folder contains data and layer files compatible with ArcGIS Pro 2.x;
  ##  - and p30 contains copies compatible with ArcGIS Pro 3.x;
  archive::archive_extract(q_map, q_map_dir)
}

pine_gdb <- file.path(q_map_dir, "p30", "pine.gdb") ## raster layer

# mdb files and extracted csvs ----------------------------------------------------------------

extracted_mdb_tables_zip <- file.path(outputPath, "extracted_mdb_tables.zip")

local({
  if (!file.exists(extracted_mdb_tables_zip)) {
    googledrive::as_id("16OK48u6g5-JdWvEFhIJ0vMxvcK-k6z5l") |>
      googledrive::drive_download(path = extracted_mdb_tables_zip, overwrite = TRUE)
    archive::archive_extract(extracted_mdb_tables_zip, outputPath)
  }
})
