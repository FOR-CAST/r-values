library(dplyr)
library(ggplot2)
library(sf)

## paths
dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
figPath <- "figures" |> fs::dir_create()
outputPath <- "outputs" |> fs::dir_create()

# load MPB SSI layers ---------------------------------------------------------

ssi_gdb <- file.path(dataPath, "MPB_SSI.gdb")

if (!file.exists(ssi_gdb)) {
  archive::archive_extract(ssi_gdb, dataPath)
}

ssi_gdb_folder <- file.path(dataPath, "MPB_SSI.gdb")
if (!dir.exists(ssi_gdb_folder)) {
  archive::archive_extract(ssi_gdb, dir = dataPath)
}

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

mpb_ssi_2008 <- get_SSI(dsn = ssi_gdb, year = 2008)
mpb_ssi_2016 <- get_SSI(dsn = ssi_gdb, year = 2016)
mpb_ssi_2023 <- get_SSI(dsn = ssi_gdb, year = 2023)

ssi_crs <- st_crs(mpb_ssi_2023)

## requires having already run 01b-import-mdb-csv.R
#all_data_sf <- file.path(outputPath, "AB", "csv", "all_mpb_site_trees_cleaned.csv") |>

##We need to use the version that has more r-values because the provincial version has too many NAs for r_value_tree
##The plot_lat_dd need to be saved as a copy for later analysis, not just removed for geometry puposes.
all_data_sf <- file.path(outputPath, "AB", "csv", "new_r_values.csv") |>
  read.csv() |>
  dplyr::filter(!is.na(plot_lon_dd) & !is.na(plot_lat_dd)) |>
  mutate(plot_lat_dd_copy = plot_lat_dd, plot_lon_dd_copy = plot_lon_dd) |>
  st_as_sf(coords = c("plot_lon_dd", "plot_lat_dd"), crs = 4326) |>
  st_make_valid() |>
  st_transform(st_crs(ssi_crs)) |>
  st_join(mpb_ssi_2008) |>
  dplyr::rename(SSI_2008 = SSI) |>
  st_join(mpb_ssi_2016) |>
  dplyr::rename(SSI_2016 = SSI) |>
  st_join(mpb_ssi_2023) |>
  dplyr::rename(SSI_2023 = SSI)

fs::dir_create(file.path(outputPath, "AB", "gdb")) |>
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

#Using our r values
ggplot(all_data_sf) +
  geom_point(aes(x = SSI_2008, y = r))

ggplot(all_data_sf) +
  geom_point(aes(x = SSI_2016, y = r))

ggplot(all_data_sf) +
  geom_point(aes(x = SSI_2023, y = r))

ggplot(all_data_sf) +
  geom_sf(aes(color = !is.na(SSI_2008))) +
  labs(title = "Trees with SSI_2008 Join", color = "Joined")

ggplot(all_data_sf) +
  geom_sf(aes(color = !is.na(SSI_2016))) +
  labs(title = "Trees with SSI_2016 Join", color = "Joined")

ggplot(all_data_sf) +
  geom_sf(aes(color = !is.na(SSI_2023))) +
  labs(title = "Trees with SSI_2023 Join", color = "Joined")

## plot the three joins in one window
library(dplyr)
library(tidyr)
library(ggplot2)

ssi_long <- all_data_sf |>
  select(geometry, SSI_2008, SSI_2016, SSI_2023) |>
  mutate(id = row_number()) |>
  pivot_longer(cols = starts_with("SSI"), names_to = "year", values_to = "ssi") |>
  mutate(joined = !is.na(ssi))

ggplot(ssi_long) +
  geom_sf(aes(color = joined)) +
  facet_wrap(~year) +
  scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "steelblue")) +
  labs(title = "SSI Join Coverage by Year", color = "Joined")

##There are lots of FALSE joins, so we will aggressively join using "nearest feature"
# Trees with missing SSI_2008
trees_missing <- all_data_sf |> filter(is.na(SSI_2008))

# Find nearest polygon index
nearest_index <- st_nearest_feature(trees_missing, mpb_ssi_2008)

# Extract nearest polygons
nearest_polygons <- mpb_ssi_2008[nearest_index, ]

# Add SSI_2008 from nearest polygon
trees_missing$SSI_2008 <- nearest_polygons$SSI

# Combine with original trees that already had SSI_2008
all_data_sf_aggressive <- all_data_sf |>
  filter(!is.na(SSI_2008)) |>
  bind_rows(trees_missing)

ggplot(all_data_sf_aggressive) +
  geom_sf(aes(color = !is.na(SSI_2008))) +
  labs(title = "Trees with aggressive SSI_2008 Join", color = "Joined")

#2016
trees_missing <- all_data_sf |> filter(is.na(SSI_2016))

# Find nearest polygon index
nearest_index <- st_nearest_feature(trees_missing, mpb_ssi_2016)

# Extract nearest polygons
nearest_polygons <- mpb_ssi_2016[nearest_index, ]

# Add SSI_2016 from nearest polygon
trees_missing$SSI_2016 <- nearest_polygons$SSI

# Combine with original trees that already had SSI_2016
all_data_sf_aggressive <- all_data_sf |>
  filter(!is.na(SSI_2016)) |>
  bind_rows(trees_missing)

ggplot(all_data_sf_aggressive) +
  geom_sf(aes(color = !is.na(SSI_2016))) +
  labs(title = "Trees with aggressive SSI_2016 Join", color = "Joined")

#2023
trees_missing <- all_data_sf |> filter(is.na(SSI_2023))

# Find nearest polygon index
nearest_index <- st_nearest_feature(trees_missing, mpb_ssi_2023)

# Extract nearest polygons
nearest_polygons <- mpb_ssi_2023[nearest_index, ]

# Add SSI_2023 from nearest polygon
trees_missing$SSI_2023 <- nearest_polygons$SSI

# Combine with original trees that already had SSI_2008
all_data_sf_aggressive <- all_data_sf |>
  filter(!is.na(SSI_2023)) |>
  bind_rows(trees_missing)

ggplot(all_data_sf_aggressive) +
  geom_sf(aes(color = !is.na(SSI_2023))) +
  labs(title = "Trees with aggressive SSI_2023 Join", color = "Joined")

#How "aggressive" was the "aggressive join procedure?" (Test using 2008)
trees_missing$distance_to_ssi_2008 <- st_distance(trees_missing, nearest_polygons, by_element = TRUE)

ggplot(trees_missing, aes(x = as.numeric(distance_to_ssi_2008) / 1000)) +
  geom_histogram(binwidth = 1, fill = "darkorange", color = "white") +
  labs(title = "Distance to Nearest SSI Polygon (2008)", x = "Distance (km)", y = "Tree Count")

#Next piece assumes explore_Q_maps has been run
library(terra)

# Extract Q-values from raster to tree points
all_data_sf$Q <- terra::extract(pine_q, vect(all_data_sf))$pine

summary(all_data_sf$Q)
hist(all_data_sf$Q, breaks = 50, col = "skyblue", main = "Distribution of Q-values")

all_data_df <- st_drop_geometry(all_data_sf)

abr.early <- all_data_df |> filter(beetle_yr <= 2015)
abr.late  <- all_data_df |> filter(beetle_yr >= 2016)

sapply(abr.early[, c("r", "beetle_yr", "plot_lat_dd_copy", "plot_lon_dd_copy", "nbr_infested", "dbh", "ht_pitch_tube")], function(x) sum(is.na(x)))
sapply(abr.late[, c("r", "beetle_yr", "plot_lat_dd_copy", "plot_lon_dd_copy", "nbr_infested", "dbh", "ht_pitch_tube")], function(x) sum(is.na(x)))

abr.lm.early <- lm(
  log10(r + 1) ~ beetle_yr + plot_lat_dd_copy + plot_lon_dd_copy + Q +
    log10(nbr_infested + 1) * dbh * ht_pitch_tube*SSI_2008,
  data = abr.early,
  na.action = na.omit
)
summary(abr.lm.early)

abr.lm.late <- lm(
  log10(r + 1) ~ beetle_yr + plot_lat_dd_copy + plot_lon_dd_copy +
    log10(nbr_infested + 1) * dbh * ht_pitch_tube,
  data = abr.late,
  na.action = na.omit
)
summary(abr.lm.late)

abr.lm.late <- lm(
  log10(r + 1) ~ beetle_yr + plot_lat_dd_copy + plot_lon_dd_copy + Q +
    log10(nbr_infested + 1) * dbh * ht_pitch_tube * SSI_2016,
  data = abr.late,
  na.action = na.omit
)
summary(abr.lm.late)

library(car)
#examine variance inflation factor associated with multi-colinearity
vif(abr.lm.late)

