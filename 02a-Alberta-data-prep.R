## append r-values ----------------------------------------------------------------------------

sum(abr$r_value_site[!is.na(abr$r_value_site)] == -999) ## 0
sum(abr$r_value_tree[!is.na(abr$r_value_tree)] == -999) ## 0

abr$r_value_site[abr$r_value_site[!is.na(abr$r_value_site)] == -999]
abr$r_value_site[abr$r_value_site[!is.na(abr$r_value_site)] == -999] <- NA

abr$r_value_site == -999

str(abr)

table(abr$nbr_infested)

abr$live <- rowSums(
  abr[, c(
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
  )],
  na.rm = TRUE
)

dev.new()
hist(log10(abr$live + 1))

abr$holes <- rowSums(
  abr[, c(
    "ns1_holes",
    "ns2_holes",
    "ss1_holes",
    "ss2_holes"
  )],
  na.rm = TRUE
)

dev.new()
hist(log10(abr$holes + 1))

## adding 1 in denominator avoids a meaningless division by zero "error",
## with a relatively small cost in basing the r-value low
abr$r <- abr$live / (abr$holes + 1)

abr_csv <- file.path(outputPath, "AB", "csv", "new_r_values.csv")

write.csv(abr, abr_csv, row.names = FALSE)

# MPB SSI layers ------------------------------------------------------------------------------

source("R/mpb_ssi.R")

ssi_gdb <- file.path(dataPath, "MPB_SSI.gdb")

if (!(file.exists(ssi_gdb) || dir.exists(ssi_gdb))) {
  archive::archive_extract(ssi_gdb, dataPath)
}

# build model data frame ----------------------------------------------------------------------

## check if data for modelling exists. If it doesn't, build it.
model_data_csv <- file.path(outputPath, "AB", "csv", "new_r_values_w_QSSIPVOL.csv")

if (!file.exists(model_data_csv) || rerun_all) {
  local({
    mpb_ssi_2008 <- get_SSI(dsn = ssi_gdb, year = 2008)
    mpb_ssi_2016 <- get_SSI(dsn = ssi_gdb, year = 2016)
    mpb_ssi_2023 <- get_SSI(dsn = ssi_gdb, year = 2023)

    ssi_crs <- sf::st_crs(mpb_ssi_2023)

    ## use the version that has more r-values because the provincial version has too many NAs for r_value_tree
    ## keep plot_lat/lon_dd for later analysis, not just removed for geometry purposes.
    all_data_sf <- read.csv(abr_csv) |>
      filter(!is.na(plot_lon_dd) & !is.na(plot_lat_dd)) |>
      mutate(plot_lat_dd_copy = plot_lat_dd, plot_lon_dd_copy = plot_lon_dd) |>
      sf::st_as_sf(coords = c("plot_lon_dd", "plot_lat_dd"), crs = 4326) |>
      sf::st_make_valid() |>
      sf::st_transform(sf::st_crs(ssi_crs)) |>
      sf::st_join(mpb_ssi_2008) |>
      rename(SSI_2008 = SSI) |>
      sf::st_join(mpb_ssi_2016) |>
      rename(SSI_2016 = SSI) |>
      sf::st_join(mpb_ssi_2023) |>
      rename(SSI_2023 = SSI)

    fs::dir_create(file.path(outputPath, "AB", "gdb")) |>
      file.path("all_mpb_site_trees_cleaned_SSI.gdb") |>
      sf::st_write(all_data_sf, dsn = _, append = FALSE)

    ## SSI correlations
    all_data_sf |>
      sf::st_drop_geometry() |>
      select(starts_with("SSI")) |>
      na.omit() |>
      stats::cor()

    ## SSI vs r_value_tree
    if (plot_all) {
      ggplot(all_data_sf) +
        geom_point(aes(x = SSI_2008, y = r_value_tree))

      ggplot(all_data_sf) +
        geom_point(aes(x = SSI_2016, y = r_value_tree))

      ggplot(all_data_sf) +
        geom_point(aes(x = SSI_2023, y = r_value_tree))
    }

    ## SSI vs our r values
    if (plot_all) {
      ggplot(all_data_sf) +
        geom_point(aes(x = SSI_2008, y = r))

      ggplot(all_data_sf) +
        geom_point(aes(x = SSI_2016, y = r))

      ggplot(all_data_sf) +
        geom_point(aes(x = SSI_2023, y = r))

      ggplot(all_data_sf) +
        geom_sf(aes(color = !is.na(SSI_2008))) +
        geom_sf(data = ab_sf, fill = NA) +
        labs(title = "Trees with SSI_2008 Join", color = "Joined")

      ggplot(all_data_sf) +
        geom_sf(aes(color = !is.na(SSI_2016))) +
        geom_sf(data = ab_sf, fill = NA) +
        labs(title = "Trees with SSI_2016 Join", color = "Joined")

      ggplot(all_data_sf) +
        geom_sf(aes(color = !is.na(SSI_2023))) +
        geom_sf(data = ab_sf, fill = NA) +
        labs(title = "Trees with SSI_2023 Join", color = "Joined")
    }

    ## plot the three joins in one window
    ssi_long <- all_data_sf |>
      select(geometry, starts_with("SSI")) |>
      mutate(id = row_number()) |>
      pivot_longer(cols = starts_with("SSI"), names_to = "year", values_to = "ssi") |>
      mutate(joined = !is.na(ssi))

    if (plot_all) {
      # dev.new()
      ggplot(ssi_long) +
        geom_sf(aes(color = joined)) +
        facet_wrap(~year) +
        geom_sf(data = ab_sf, fill = NA) +
        scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "steelblue")) +
        labs(title = "SSI Join Coverage by Year", color = "Joined")
    }

    ## There are lots of FALSE joins, so we will aggressively join using "nearest feature"

    ## 2008 ---------------------------------------------------------------------
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

    if (plot_all) {
      ggplot(all_data_sf_aggressive) +
        geom_sf(aes(color = !is.na(SSI_2008))) +
        geom_sf(data = ab_sf, fill = NA) +
        labs(title = "Trees with aggressive SSI_2008 Join", color = "Joined")
    }

    ## 2016 ---------------------------------------------------------------------
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

    if (plot_all) {
      ggplot(all_data_sf_aggressive) +
        geom_sf(aes(color = !is.na(SSI_2016))) +
        geom_sf(data = ab_sf, fill = NA) +
        labs(title = "Trees with aggressive SSI_2016 Join", color = "Joined")
    }

    ## 2023 ---------------------------------------------------------------------
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

    if (plot_all) {
      ggplot(all_data_sf_aggressive) +
        geom_sf(aes(color = !is.na(SSI_2023))) +
        geom_sf(data = ab_sf, fill = NA) +
        labs(title = "Trees with aggressive SSI_2023 Join", color = "Joined")
    }

    ## How "aggressive" was the "aggressive join procedure?" (Test using 2008)
    trees_missing$distance_to_ssi_2008 <- st_distance(
      trees_missing,
      nearest_polygons,
      by_element = TRUE
    )

    if (plot_all) {
      # dev.new()
      ggplot(trees_missing, aes(x = as.numeric(distance_to_ssi_2008) / 1000)) +
        geom_histogram(binwidth = 1, fill = "darkorange", color = "white") +
        labs(
          title = "Distance to Nearest SSI Polygon (2008)",
          x = "Distance (km)",
          y = "Tree Count"
        )
    }
    ## Extract Q-values from raster to tree points
    all_data_sf$Q <- terra::extract(pine_q, terra::vect(all_data_sf))$pine

    summary(all_data_sf$Q)

    if (plot_all) {
      # dev.new()
      hist(all_data_sf$Q, breaks = 50, col = "skyblue", main = "Distribution of Q-values")
    }

    ## Extract pine volume per hectare
    all_data_sf$PineVol <- terra::extract(
      x = bleiker2019,
      y = terra::vect(all_data_sf)
    )$Overstory_Raster_PineVol

    summary(all_data_sf$PineVol)

    ## Prepare data for modelling
    all_data_df <- st_drop_geometry(all_data_sf)

    ## Save to disk
    write.csv(all_data_df, model_data_csv, row.names = FALSE)
  })
}

all_data_df <- read.csv(model_data_csv)

# get MPB winter mortality (winterkill) -------------------------------------------------------

source("R/biosim.R")

if (FALSE) {
  BioSIM::getModelList() ## list the models available

  BioSIM::getModelHelp("MPB_Cold_Tolerance_Annual")
}

all_data_rds <- file.path(outputPath, "AB", "all_data_df_clean.rds")

if (!file.exists(all_data_rds) || rerun_all) {
  local({
    ## generate an elevation for every location in all_data_df
    coords <- data.frame(x = all_data_df$plot_lon_dd_copy, y = all_data_df$plot_lat_dd_copy)
    names(coords) <- c("lon", "lat")
    coords_sf <- st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
    elevations <- elevatr::get_elev_point(
      locations = coords_sf,
      prj = "+proj=longlat +datum=WGS84",
      src = "aws"
    )

    all_data_df$elevation <- elevations$elevation

    ## rename plot_lat/lon_dd_copy to simply lat/lon
    all_data_df <- all_data_df |>
      rename(
        lon = plot_lon_dd_copy,
        lat = plot_lat_dd_copy
      )

    saveRDS(all_data_df, all_data_rds)
  })
}

all_data_df <- readRDS(all_data_rds)

## test on just the unique locations; generate a unique list for BioSIM

site_year_df <- all_data_df |>
  select(lat, lon, beetle_yr, elevation) |>
  distinct()

site_year__MPBwk_rds <- file.path(outputPath, "AB", "site_year__MPBwk_results.rds")

if (file.exists(site_year__MPBwk_rds) && !rerun_all) {
  site_year__MPBwk_results <- readRDS(site_year__MPBwk_rds)
} else {
  site_year__MPBwk_results <- mpb_cold_tol(site_year_df)
  site_year__MPBwk_results <- site_year__MPBwk_results |>
    mutate(Psurv_prop = Psurv / 100)

  saveRDS(site_year__MPBwk_results, site_year__MPBwk_rds)
}

plot(site_year__MPBwk_results$Tmin, site_year__MPBwk_results$Psurv)

## verify results by mapping
site_year_sf <- st_as_sf(
  site_year__MPBwk_results,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)
ggplot(site_year_sf) +
  geom_sf(aes(color = Psurv), size = 2) +
  geom_sf(data = ab_sf, fill = NA) +
  scale_color_viridis_c(option = "plasma", name = "Survival (%)") +
  facet_wrap(~Year) +
  theme_minimal() +
  labs(
    title = "Winter Survival Probability by Site and Year",
    subtitle = "BioSIM MPB Cold Tolerance Model",
    caption = "Each point represents a unique site-year combination"
  )

psurv_summary <- site_year__MPBwk_results |>
  group_by(Year) |>
  summarise(
    mean_Psurv = mean(Psurv),
    sd_Psurv = sd(Psurv),
    n = n()
  )

tmin_summary <- site_year__MPBwk_results |>
  group_by(Year) |>
  summarise(
    mean_Tmin = mean(Tmin),
    sd_Tmin = sd(Tmin),
    n = n()
  )

## plot Psurv and Tmin over time

gg_psurv_summary <- ggplot(psurv_summary, aes(x = Year, y = mean_Psurv)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue", size = 2) +
  geom_ribbon(
    aes(ymin = mean_Psurv - sd_Psurv, ymax = mean_Psurv + sd_Psurv),
    alpha = 0.2,
    fill = "blue"
  ) +
  labs(
    title = "Mean Winter Survival Probability Over Time",
    y = "Mean Psurv (%)",
    x = "Year",
    caption = "Shaded area shows ±1 SD across sites"
  ) +
  theme_minimal()

ggsave(file.path(figPath, "mean_Psurv_over_time.png"), gg_psurv_summary)

gg_tmin_summary <- ggplot(tmin_summary, aes(x = Year, y = mean_Tmin)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue", size = 2) +
  geom_ribbon(
    aes(ymin = mean_Tmin - sd_Tmin, ymax = mean_Tmin + sd_Tmin),
    alpha = 0.2,
    fill = "blue"
  ) +
  labs(
    title = "Mean Minimum Winter Temperature Over Time",
    y = "Mean Tmin",
    x = "Year",
    caption = "Shaded area shows ±1 SD across sites"
  ) +
  theme_minimal()

ggsave(file.path(figPath, "mean_Tmin_over_time.png"), gg_tmin_summary)

## run on all 13312 samples

all_data_df_join_Psurv_csv <- file.path(outputPath, "AB", "csv", "new_r_values_w_Q_SSI_P_PVOL.csv")

if (file.exists(all_data_df_join_Psurv_csv) && !rerun_all) {
  all_data_df_join_Psurv <- read.csv(all_data_df_join_Psurv_csv)
} else {
  site_year_results <- mpb_cold_tol(all_data_df)

  site_year_results_min <- site_year_results |>
    select(row_index, Psurv, Tmin)

  all_data_df_join_Psurv <- all_data_df |>
    mutate(row_index = row_number()) |>
    left_join(site_year_results_min, by = "row_index")

  str(all_data_df_join_Psurv)

  write.csv(all_data_df_join_Psurv, all_data_df_join_Psurv_csv, row.names = FALSE)
}

# get Climate Moisture Index (CMI) ------------------------------------------------------------

if (FALSE) {
  BioSIM::getModelList() ## list the models available

  BioSIM::getModelHelp("Climate_Mosture_Index_Annual")
}

all_data_df_join_CMI_csv <- file.path(
  outputPath,
  "AB",
  "csv",
  "new_r_values_w_Q_SSI_P_PVOL_CMI.csv"
)

if (!file.exists(all_data_df_join_CMI_csv) || rerun_all) {
  local({
    site_year_results <- biosim_cmi(all_data_df)

    site_year_results_min <- site_year_results |>
      select(row_index, CMI)

    all_data_df_join_CMI <- all_data_df_join_Psurv |>
      mutate(row_index = row_number()) |>
      left_join(site_year_results_min, by = "row_index")

    str(all_data_df_join_CMI)

    write.csv(all_data_df_join_CMI, all_data_df_join_CMI_csv, row.names = FALSE)
  })
}

all_data_df_join_CMI <- read.csv(all_data_df_join_CMI_csv)
