# packages ------------------------------------------------------------------------------------

library(DBI)
library(odbc)
library(ggplot2)
library(ggspatial)

# check if data for modelling exists. If it doesn't, build it. ---------------------------------------------------------

model.data <- file.path(outputPath, "AB", "csv", "new_r_values_w_QSSI.csv")

if (!file.exists(model.data)) {

  # read in data from MS Access databases -------------------------------------------------------

  ## TODO: see 01b-import-mdb-csv.R

  # setup ---------------------------------------------------------------------------------------

  ## paths
  dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
  figPath <- "figures" |> fs::dir_create()
  outputPath <- "outputs" |> fs::dir_create()

  #Read in site/tree data to which we will append our calculations of r-value
  abr <- read.csv(file.path(outputPath, "AB", "csv", "all_mpb_site_trees_cleaned.csv"), header = TRUE)

  ## data check
  dev.new()
  par(mfrow = c(2, 4))
  plot(abr$lon_dd, abr$lat_dd, col = "red")
  plot(abr$plot_lon_dd, abr$plot_lat_dd, col = "blue")
  hist(abr$beetle_yr, breaks = c(2006:2019))
  hist(log10(abr$nbr_infested + 1))
  hist(abr$dbh[abr$dbh < 100])
  hist(abr$ht_pitch_tube)
  hist(log10(abr$r_value_site) + 1)
  hist(log10(abr$r_value_tree) + 1)

  h <- hist(abr$nbr_infested, plot = FALSE)
  plot(h, yaxt = "n", main = "Histogram with Log-Scaled Y-Axis", xlab = "nbr_infested")
  axis(2, at = pretty(log10(h$counts)), labels = 10^pretty(log10(h$counts)))

  table(abr$beetle_yr)

  hist(abr$dbh[!is.na(abr$dbh)])

  sum(abr$r_value_site[!is.na(abr$r_value_site)] == -999)
  sum(abr$r_value_tree[!is.na(abr$r_value_tree)] == -999)

  abr$r_value_site[all_data$r_value_site[!is.na(all_data$r_value_site)] == -999]

  abr$r_value_site[all_data$r_value_site[!is.na(all_data$r_value_site)] == -999] <- NA

  abr$r_value_site == -999

  str(abr)

  table(abr$nbr_infested)

  library(dplyr)

  abr$live <- rowSums(abr[, c(
    "ns1_eggs_live",
    "ns1_larvae_live",
    "ns1_pupae_live",
    "ns1_teneral_adults_live",
    "ns2_eggs_live",
    "ns2_larvae_live",
    "ns2_pupae_live",
    "ns2_teneral_adults_live",
    "ss1_eggs_live",
    "ss1_larvae_live",
    "ss1_pupae_live",
    "ss1_teneral_adults_live",
    "ss2_eggs_live",
    "ss2_larvae_live",
    "ss2_pupae_live",
    "ss2_teneral_adults_live"
  )], na.rm = TRUE)

  dev.new()
  hist(log10(abr$live + 1))

  abr$holes <- rowSums(abr[, c(
    "ns1_holes",
    "ns2_holes",
    "ss1_holes",
    "ss2_holes"
  )], na.rm = TRUE)
  dev.new()
  hist(log10(abr$holes + 1))

  ## adding 1 in denominator avoids a meaningless division by zero "error",
  ## with a relatively small cost in basing the r-value low
  abr$r <- abr$live / (abr$holes + 1)

  ## Write new dataset containing the r values
  write.csv(abr, file.path(outputPath, "AB", "csv", "new_r_values.csv"), row.names = FALSE)

  ############
  # SSI      #
  ############

  #Borrowed from script explore_SSI

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

  ## We need to use the version that has more r-values because the provincial version has too many NAs for r_value_tree
  ## The plot_lat_dd need to be saved as a copy for later analysis, not just removed for geometry purposes.
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

  dev.new()
  ggplot(ssi_long) +
    geom_sf(aes(color = joined)) +
    facet_wrap(~year) +
    scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "steelblue")) +
    labs(title = "SSI Join Coverage by Year", color = "Joined")

  ##There are lots of FALSE joins, so we will aggressively join using "nearest feature"
  ## Trees with missing SSI_2008
  trees_missing <- all_data_sf |> filter(is.na(SSI_2008))

  ## Find nearest polygon index
  nearest_index <- st_nearest_feature(trees_missing, mpb_ssi_2008)

  ## Extract nearest polygons
  nearest_polygons <- mpb_ssi_2008[nearest_index, ]

  ## Add SSI_2008 from nearest polygon
  trees_missing$SSI_2008 <- nearest_polygons$SSI

  ## Combine with original trees that already had SSI_2008
  all_data_sf_aggressive <- all_data_sf |>
    filter(!is.na(SSI_2008)) |>
    bind_rows(trees_missing)

  ggplot(all_data_sf_aggressive) +
    geom_sf(aes(color = !is.na(SSI_2008))) +
    labs(title = "Trees with aggressive SSI_2008 Join", color = "Joined")

  #2016
  trees_missing <- all_data_sf |> filter(is.na(SSI_2016))

  ## Find nearest polygon index
  nearest_index <- st_nearest_feature(trees_missing, mpb_ssi_2016)

  ## Extract nearest polygons
  nearest_polygons <- mpb_ssi_2016[nearest_index, ]

  ## Add SSI_2016 from nearest polygon
  trees_missing$SSI_2016 <- nearest_polygons$SSI

  ## Combine with original trees that already had SSI_2016
  all_data_sf_aggressive <- all_data_sf |>
    filter(!is.na(SSI_2016)) |>
    bind_rows(trees_missing)

  ggplot(all_data_sf_aggressive) +
    geom_sf(aes(color = !is.na(SSI_2016))) +
    labs(title = "Trees with aggressive SSI_2016 Join", color = "Joined")

  #2023
  trees_missing <- all_data_sf |> filter(is.na(SSI_2023))

  ## Find nearest polygon index
  nearest_index <- st_nearest_feature(trees_missing, mpb_ssi_2023)

  ## Extract nearest polygons
  nearest_polygons <- mpb_ssi_2023[nearest_index, ]

  ## Add SSI_2023 from nearest polygon
  trees_missing$SSI_2023 <- nearest_polygons$SSI

  ## Combine with original trees that already had SSI_2008
  all_data_sf_aggressive <- all_data_sf |>
    filter(!is.na(SSI_2023)) |>
    bind_rows(trees_missing)

  ggplot(all_data_sf_aggressive) +
    geom_sf(aes(color = !is.na(SSI_2023))) +
    labs(title = "Trees with aggressive SSI_2023 Join", color = "Joined")

  ## How "aggressive" was the "aggressive join procedure?" (Test using 2008)
  trees_missing$distance_to_ssi_2008 <- st_distance(trees_missing, nearest_polygons, by_element = TRUE)

  dev.new()
  ggplot(trees_missing, aes(x = as.numeric(distance_to_ssi_2008) / 1000)) +
    geom_histogram(binwidth = 1, fill = "darkorange", color = "white") +
    labs(title = "Distance to Nearest SSI Polygon (2008)", x = "Distance (km)", y = "Tree Count")

  ############
  # Q values #
  ############

  #Borrowed from script explore_Q_maps

  library(terra)

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
  dev.new()
  terra::plot(pine_q)

  # Extract Q-values from raster to tree points

  all_data_sf$Q <- terra::extract(pine_q, vect(all_data_sf))$pine

  summary(all_data_sf$Q)
  dev.new()
  hist(all_data_sf$Q, breaks = 50, col = "skyblue", main = "Distribution of Q-values")

  ## Prepare data for modeling

  all_data_df <- st_drop_geometry(all_data_sf)

  ## Save to disk
  write.csv(all_data_df, file.path(outputPath, "AB", "csv", "new_r_values_w_QSSI.csv"), row.names = FALSE)

} #end of data creation step

# Data exists, either in memory or on disk, so read it in.
all_data_df <- read.csv(model.data)

## switch terminology from "all_data" to "abr"
## split dataset based on 2015-2016 as pivot
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
    log10(nbr_infested + 1) * dbh * ht_pitch_tube * SSI_2016,
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

## Examine variance inflation factor associated with multi-colinearity
library(car)
vif(abr.lm.early)
vif(abr.lm.late)

## We are drowning in multi-colinearity due to SSI and Q and lat/lon.
## Try dropping lat/lon

abr.lm.early <- lm(
  log10(r + 1) ~ beetle_yr + Q +
    log10(nbr_infested + 1) * dbh * ht_pitch_tube * SSI_2008,
  data = abr.early,
  na.action = na.omit
)
summary(abr.lm.early)

abr.lm.late <- lm(
  log10(r + 1) ~ beetle_yr + Q +
    log10(nbr_infested + 1) * dbh * ht_pitch_tube * SSI_2016,
  data = abr.late,
  na.action = na.omit
)
summary(abr.lm.late)

## We are still drowning in multi-colinearity due to SSI and Q.
## Try dropping Q.
abr.lm.early <- lm(
  log10(r + 1) ~ beetle_yr +
    log10(nbr_infested + 1) * dbh * ht_pitch_tube * SSI_2008,
  data = abr.early,
  na.action = na.omit
)
summary(abr.lm.early)

abr.lm.late <- lm(
  log10(r + 1) ~ beetle_yr +
    log10(nbr_infested + 1) * dbh * ht_pitch_tube * SSI_2016,
  data = abr.late,
  na.action = na.omit
)
summary(abr.lm.late)

## We are still drowning in multi-colinearity due to 4-way interaction.
## Try dropping SSI.
abr.lm.early <- lm(
  log10(r + 1) ~ beetle_yr +
    log10(nbr_infested + 1) * dbh * ht_pitch_tube,
  data = abr.early,
  na.action = na.omit
)
summary(abr.lm.early)

abr.lm.late <- lm(
  log10(r + 1) ~ beetle_yr +
    log10(nbr_infested + 1) * dbh * ht_pitch_tube,
  data = abr.late,
  na.action = na.omit
)
summary(abr.lm.late)

## What if SSI affects immigration behaviour, not beetle pressure/colonization success?
## Try including SSI not as interaction.
abr.lm.early <- lm(
  log10(r + 1) ~ beetle_yr + SSI_2008 +
    log10(nbr_infested + 1) * dbh * ht_pitch_tube,
  data = abr.early,
  na.action = na.omit
)
summary(abr.lm.early)

abr.lm.late <- lm(
  log10(r + 1) ~ beetle_yr + SSI_2016 +
    log10(nbr_infested + 1) * dbh * ht_pitch_tube,
  data = abr.late,
  na.action = na.omit
)
summary(abr.lm.late)

## SSI ns as main effect. So remove completely and try Q.
abr.lm.early <- lm(
  log10(r + 1) ~ beetle_yr + Q +
    log10(nbr_infested + 1) * dbh * ht_pitch_tube,
  data = abr.early,
  na.action = na.omit
)
summary(abr.lm.early)

abr.lm.late <- lm(
  log10(r + 1) ~ beetle_yr + Q +
    log10(nbr_infested + 1) * dbh * ht_pitch_tube,
  data = abr.late,
  na.action = na.omit
)
summary(abr.lm.late)

vif(abr.lm.early)
vif(abr.lm.late)

## Even with the simple model the VIFs are very high, especially on HT and DBH. A shame because the DBH and HT distributions are nice.
dev.new()
HT.hist<-hist(all_data_df$ht_pitch_tube)
dev.new()
DBH.hist<-hist(all_data_df$dbh)

library(mgcv)
gam_model.e <- gam(r ~ s(dbh) + s(ht_pitch_tube) + s(log10(nbr_infested + 1)) + Q + beetle_yr, data = abr.early)
summary(gam_model.e)
gam_model.l <- gam(r ~ s(dbh) + s(ht_pitch_tube) + s(log10(nbr_infested + 1)) + Q + beetle_yr, data = abr.late)
summary(gam_model.l)

## End of testing initial models
## bring in WK (winterkill) to be computed via BioSIM API


#----------

## Setup (only need to run once on a new machine)
Sys.setenv(R_LIBCURL_SSL_REVOKE_BEST_EFFORT=TRUE) #Needed for "several things that we do in the NRCan network"

# Install j4r
install.packages("https://sourceforge.net/projects/repiceasource/files/latest/download", repos = NULL,  type="source")

# Note: On NRCan machines, must install Perforce OpenLogic OpnJDK 64-bit, as a replacement for Java
Sys.setenv(PATH = paste("C:/Program Files/OpenLogic/jdk-22.0.2.9-hotspot/bin", Sys.getenv("PATH"), sep = ";"))
options(java.home = "C:/Program Files/OpenLogic/jdk-22.0.2.9-hotspot")

# Install BioSIM
install.packages("https://sourceforge.net/projects/biosimclient.mrnfforesttools.p/files/latest/download", repos = NULL,  type="source")

#To get past NRCan firewall:
install.packages("remotes")
remotes::install_github("RNCan/BioSimClient_R")

# Load libraries:
library(BioSIM)

#--------------

# To list the models available:
BioSIM::getModelList()

