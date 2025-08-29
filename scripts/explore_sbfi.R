## NTEMS SBFI (Matasci et al. 2018)

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

url_sbfi <- "https://opendata.nfis.org/downloads/forest_change/CA_Forest_Satellite_Based_Inventory_2020.zip"
zip_sbfi <- file.path(dataPath, basename(url_sbfi))

if (!file.exists(zip_sbfi)) {
  download.file(url_sbfi, zip_sbfi)
  archive::archive_extract(zip_sbfi, dataPath)
}

# st_layers(gdb_sbfi2020)

sbfi_grid <- st_read(gdb_sbfi2020, layer = "_Grid_forested_ecosystem")

sbfi_grid_ab <- st_intersection(sbfi_grid, st_transform(ab_sf, crs(sbfi_grid)))

plot(sbfi_grid_ab["Id"])

sbfi_layer_names_ab <- paste0("CA_2020_tile_", sbfi_grid_ab[["Id"]], "_metrics")

# lyr <- sbfi_layer_names_ab[1]
lapply(sbfi_layer_names_ab, function(lyr) {
  lyr_sf <- st_read(gdb_sbfi2020, layer = lyr) |>
    select(starts_with(paste0("SPECIES_", 1:5)), starts_with("STRUCTURE_VOLUME_")) |>
    filter(
      SPECIES_1 == "PINU.CON" |
        SPECIES_1 == "PINU.BAN" |
        SPECIES_2 == "PINU.CON" |
        SPECIES_2 == "PINU.BAN" |
        SPECIES_3 == "PINU.CON" |
        SPECIES_3 == "PINU.BAN" |
        SPECIES_4 == "PINU.CON" |
        SPECIES_4 == "PINU.BAN" |
        SPECIES_5 == "PINU.CON" |
        SPECIES_5 == "PINU.BAN"
    )

  ## TODO: apportion STRUCTURE_VOLUME_TOTAL based on SPECIES_{N}_PERC for LP and JP
})

## TODO: extract pine vols for specific points for use with GAM
