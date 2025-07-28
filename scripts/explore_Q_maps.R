## introgression data from:
##   Burns, I., James, P. M., Coltman, D. W., & Cullingham, C. I. (2019).
##   Spatial and genetic structure of the lodgepole × jack pine hybrid zone.
##   Canadian Journal of Forest Research, 49(7), 844-853.

# library(archive)
# library(fs)
# library(googledrive)
# library(sf)
# limrary(terra)

## paths
dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
figPath <- "figures" |> fs::dir_create()
outputPath <- "outputs" |> fs::dir_create()

## get the data

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

pine_gdb <- file.path(q_map_dir, "p30", "pine.gdb") ## raster layer, so use terra
pine_q <- terra::rast(pine_gdb)
terra::plot(pine_q)
