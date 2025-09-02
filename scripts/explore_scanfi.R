##  Dataset:
##  Guindon L., Villemaire P., Correia D.L.P., Manka F., Lacarte S., Smiley B. 2023.
##    SCANFI: Spatialized Canadian National Forest Inventory data product. Natural Resources Canada, Canadian Forest Service, Laurentian Forestry Centre, Quebec, Canada.
##    https://doi.org/10.23687/18e6a919-53fd-41ce-b4e2-44a9707c52dc
##
##  Scientific publication:
##  Guindon L., Manka F, Correia L.P. D., Villemaire P., Smiley B., Bernier P., Gauthier S.,
##    Beaudoin A., Boucher J., Boulanger Y. A new approach for Spatializing the Canadian National
##    Forest Inventory (SCANFI) using Landsat dense time series.
##    Canadian Journal of Forest Research 2024. [In Press]

library(dplyr)
library(sf)
library(tidyr)

dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
figPath <- "figures" |> fs::dir_create()
outputPath <- "outputs" |> fs::dir_create()

ab_sf <- try(
  geodata::gadm("CAN", level = 1, path = dataPath) |>
    sf::st_as_sf() |>
    filter(NAME_1 == "Alberta") |>
    sf::st_geometry()
)

url_scanfi <- "https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/SCANFI/v1/"
files_scanfi <- file.path(
  dataPath,
  c(
    "SCANFI_sps_lodgepolePine_SW_2020_v1.2.tif",
    "SCANFI_sps_lodgepolePine_SW_2020_v1.2.tif.aux.xml",
    "SCANFI_sps_lodgepolePine_SW_2020_v1.2.tif.ovr",
    "SCANFI_sps_jackPine_SW_2020_v1.2.tif",
    "SCANFI_sps_jackPine_SW_2020_v1.2.tif.aux.xml",
    "SCANFI_sps_jackPine_SW_2020_v1.2.tif.ovr"
  )
)

if (!all(file.exists(files_scanfi))) {
  download.file(paste0(url_scanfi, basename(files_scanfi)), files_scanfi)
}

## percent cover pine (not volume)
scanfi <- file.path(files_scanfi) |>
  grep("[.]tif$", x = _, value = TRUE) |>
  terra::rast()
scanfi <- terra::crop(
  x = scanfi,
  y = sf::st_transform(ab_sf, terra::crs(scanfi))
)
