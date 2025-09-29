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

## National Parks ------------------------------------------------------------------------------

latlon <- crs("epsg:4326")
targetCRS <- crs(paste(
  "+proj=aea +lat_1=49 +lat_2=67 +lat_0=0 +lon_0=-112 +x_0=0 +y_0=0",
  "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
))

## Banff and Jasper National Parks
## see: https://hub.arcgis.com/datasets/dd8cd91871534c9aa34310eed84fe076_1/about
np_url <- "https://drive.google.com/file/d/1Rz8BzyWtirXuCbAxx8MqeJs2KlNYO9eh/"
np_file <- "National_Parks_and_National_Park_Reserves_of_Canada_Legislative_Boundaries"
np_fext <- c("cpg", "dbf", "prj", "shp", "shx")
np_zip <- file.path(dataPath, paste0(np_file, ".zip"))

if (!file.exists(np_zip)) {
  googledrive::drive_download(googledrive::as_id(np_url), np_zip)
}

if (!all(file.exists(file.path(dataPath, paste0(np_file, ".", np_fext))))) {
  archive::archive_extract(np_zip, dataPath)
}

natl_prks.latlon <- st_read(file.path(dataPath, paste0(np_file, ".shp")))
natl_prks <- st_transform(natl_prks.latlon, targetCRS)

np_banff.latlon <- natl_prks.latlon[
  natl_prks.latlon$adminAreaN == "BANFF NATIONAL PARK OF CANADA",
]
np_banff <- natl_prks[natl_prks$adminAreaN == "BANFF NATIONAL PARK OF CANADA", ]

np_jasper.latlon <- natl_prks.latlon[
  natl_prks.latlon$adminAreaN == "JASPER NATIONAL PARK OF CANADA",
]
np_jasper <- natl_prks[natl_prks$adminAreaN == "JASPER NATIONAL PARK OF CANADA", ]

## Combine parks
parks <- rbind(np_banff, np_jasper)

## Get bounding box and convert to sf polygon
bbox_parks <- st_bbox(parks) |> st_as_sfc() |> st_sf()

## highways and roads -------------------------------------------------------------------------

nrn_dir <- file.path(dataPath, "NRCan") |> fs::dir_create()
nrn_url <- "https://geo.statcan.gc.ca/nrn_rrn/ab/nrn_rrn_ab_SHAPE.zip"
nrn_zip <- file.path(dirname(nrn_dir), basename(nrn_url))

if (!file.exists(nrn_zip)) {
  download.file(nrn_url, destfile = nrn_zip)
  archive::archive_extract(nrn_zip, nrn_dir)
}

ab_roads <- file.path(
  nrn_dir,
  "NRN_RRN_AB_SHAPE",
  "NRN_AB_17_0_SHAPE_en",
  "NRN_AB_17_0_ROADSEG.shp"
) |>
  st_read() |>
  st_transform(st_crs(parks))

## Filter for major highways only
ab_roads_clipped <- ab_roads |>
  filter(
    RTNUMBER1 %in%
      c("1", "16", "93") |
      grepl("Highway 1|Highway 16|Highway 93", R_STNAME_C, ignore.case = TRUE)
  ) |>
  st_intersection(bbox_parks)

## water bodies -------------------------------------------------------------------------------

## see https://ftp.maps.canada.ca/pub/nrcan_rncan/vector/canvec/fgdb/Hydro/
ab_hydro_dir <- file.path(dataPath, "NRCan", "canvec_250K_AB_Hydro") |> fs::dir_create()
ab_hydro_url <- "https://ftp.maps.canada.ca/pub/nrcan_rncan/vector/canvec/fgdb/Hydro/canvec_250K_AB_Hydro_fgdb.zip"

ab_hydro_zip <- file.path(dirname(ab_hydro_dir), basename(ab_hydro_url))
ab_hydro_gdb <- file.path(ab_hydro_dir, "canvec_250K_AB_Hydro.gdb")

if (!file.exists(ab_hydro_zip)) {
  download.file(ab_hydro_url, destfile = ab_hydro_zip)
}

if (!file.exists(ab_hydro_gdb)) {
  archive::archive_extract(ab_hydro_zip, dir = ab_hydro_dir)
}

ab_hydro <- sf::st_read(ab_hydro_gdb, layer = "waterbody_2") |>
  sf::st_transform(st_crs(parks))

ab_hydro_clipped <- st_intersection(ab_hydro, bbox_parks)

# pine maps -----------------------------------------------------------------------------------

gdb_bleiker2019 <- file.path(dataPath, "AB_PineVolumes_Lambert.gdb")

local({
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

# MPB SSI layers ------------------------------------------------------------------------------

source("R/mpb_ssi.R") ## helper fun to load the layers

ssi_gdb <- file.path(dataPath, "MPB_SSI.gdb")
ssi_zip <- paste0(ssi_gdb, ".zip")

if (!(file.exists(ssi_gdb) || dir.exists(ssi_gdb))) {
  googledrive::as_id("1meuraFVblxD8ZlwuqgKy0ADHOs-EAN13") |>
    googledrive::drive_download(path = ssi_zip)
  archive::archive_extract(ssi_zip, dataPath)
}

# Olesinski Jasper Banff 2012-2023 ------------------------------------------------------------

MPB_Jasper_Banff_2012_2023_zip <- file.path(dataPath, "MPB_Jasper_Banff_2012_2023.zip")
MPB_Jasper_Banff_2012_2023_shp <- file.path(
  dataPath,
  "MPB_Jasper_Banff_2012_2023",
  "MPB_Jasper_Banff_2012_2023.shp"
)

local({
  if (!file.exists(MPB_Jasper_Banff_2012_2023_zip)) {
    googledrive::as_id("1TwBJyv-7lgR2GJEN098yYQuflRO5hw_D") |>
      googledrive::drive_download(path = MPB_Jasper_Banff_2012_2023_zip)
    archive::archive_extract(MPB_Jasper_Banff_2012_2023_zip, dataPath)
  }
})
