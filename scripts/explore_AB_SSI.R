library(dplyr)
library(ggplot2)
library(sf)

## paths
dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
figPath <- "figures" |> fs::dir_create()
outputPath <- "outputs" |> fs::dir_create()

## load MPB SSI layers
get_SSI <- function(dsn, year) {
  ssi <- st_read(dsn = dsn, layer = paste0("MPB_SSI_", year)) |>
    st_cast("MULTIPOLYGON")

  ssi <- st_make_valid(ssi)

  ## keep only the SSI values and polygon geometries,
  ## and use 'SSI' as the column/field name
  ssi <- ssi |>
    select(matches("^(MPB_SSI|SSI)$"), Shape) |>
    rename(any_of(c(SSI = "MPB_SSI")))

  return(ssi)
}

mpb_ssi_2008 <- get_SSI(dsn = "data/MPB_SSI.gdb", year = 2008)
mpb_ssi_2016 <- get_SSI(dsn = "data/MPB_SSI.gdb", year = 2016)
mpb_ssi_2023 <- get_SSI(dsn = "data/MPB_SSI.gdb", year = 2023)

ssi_crs <- st_crs(mpb_ssi_2023)

## requires having already run 01b-import-mdb-csv.R
all_data_sf <- file.path(dataPath, "AB", "csv", "all_mpb_site_trees_cleaned.csv") |>
  read.csv() |>
  dplyr::filter(!is.na(lon_dd) & !is.na(lat_dd)) |>
  st_as_sf(coords = c("lon_dd","lat_dd"), crs = 4326) |>
  st_make_valid() |>
  st_transform(st_crs(ssi_crs)) |>
  st_join(mpb_ssi_2008) |>
  dplyr::rename(SSI_2008 = SSI) |>
  st_join(mpb_ssi_2016) |>
  dplyr::rename(SSI_2016 = SSI) |>
  st_join(mpb_ssi_2023) |>
  dplyr::rename(SSI_2023 = SSI)

fs::dir_create(file.path(dataPath, "AB", "gdb")) |>
  file.path("all_mpb_site_trees_cleaned_SSI.gdb") |>
  st_write(all_data_sf, dsn = _)

## SSI correlations

all_data_sf |>
  st_drop_geometry() |>
  dplyr::select(starts_with("SSI")) |>
  na.omit() |>
  cor()

## SSI vs r_value_tree

ggplot(all_data_sf) +
  geom_point(aes(x = SSI_2008, y = r_value_tree))

ggplot(all_data_sf) +
  geom_point(aes(x = SSI_2016, y = r_value_tree))

ggplot(all_data_sf) +
  geom_point(aes(x = SSI_2023, y = r_value_tree))
