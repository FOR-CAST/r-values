# In this script we examine three components of mountain pine beetle population
# dynamics in two Rocky mountain parks, Jasper and Banff, where no control was undertaken.
# There are more data for Jasper than Banff because the outbreak was more intense
# at Banff, so the analysis is asymmetric. The three data components are:
# 1. trends in red tree counts/area (capital "R" is the interannual rate of change in Xt+1/Xt);
# 2. r-values (brood productivity) (lowercase "r");
# 3. predicted winter mortality (1-Psurv) (and drought; CMI) (using BioSIM).
#
# Our primary question is whether (3) predicts (2) predicts (1).
#
# 1. Tree counts. For both Jasper and Banff we have red tree areas 2013-2021, the peak outbreak years.
# Prior to outbreak, during 1999-2012, we have tree counts for Banff. We have them for Jasper too,
# but 2011 and 2012 are missing, although there is a count for 2013.
#
# 2. r-values. These are based on 4" disks that average just 0.5 female entrance holes per disk.
# There are no r-values for Banff. For Jasper there are r-values for beetle years
# 2014, 2015, 2016, recorded by ASRD in .mdb files. The survey is done in the year after the "beetle year"
# These need to be extracted, converted to .csv tables, and site and tree data merged
# (as we did for the rest of Alberta in a companion paper). The formats of the three .mdb files are the same
# (and very similar to the provincial files for 2006-2019 analyzed in the companion paper).
# The r-values for beetle years 2017, 2018, 2020, 2021 are in .shp files. There are no data for 2019 because of covid
# during the sampling year, 2020. These files are not named for beetle year, but survey year.
# The format of the .shp attribute tables varies somewhat.
#
# 3. We run BioSIM winter mortality predictions for all the years and locations for which we have comparable data.
# The output variable is called Psurv - the probability of surviving the winter. The model takes initial phenology
# of life stages as initial conditions, but we always use the default assumptions of a stable seasonality, because
# our population samples are so small (based on 4" disks) we can't reliably infer phenology.
# There are two sets of simulations to be run. The first is for the townsites of Jasper and Banff, which have weather stations,
# for the beetle years 1998-2022, which fully brackets the years of red tree counts and r-values. We expect this to show a
# long warming trend in winter weather, followed by successive cold snaps in the last three years, when the outbreak was dying off.
# The second is for the specific locations in JNP for which we have R-values, 2014-2021. These simulations will be used to create
# data tables for modelling, and yearly maps for the whole Jasper NP, for presentation purposes. Those maps will likely be
# aggregated into a single pair of maps representing the average during the rise, and the average during the collapse. BioSIM does
# spatial co-kriging efficiently, but we can't use that functionality in the web API. So we work around that here, in R,
# by sending a high density grid of points across JNP and BNP to simulate. Maps are smoothed over the lat and lon grid.
# We also run the Climate Moisture Index (CMI) model in BioSIM, as dryness has been fingered as a key determinant of outbreak potential,
# and just as winter temperature is known to fluctuate severely across years, so does drought.
#
# For the beetle years 2014-16 the Jasper r-values data are "rich" (as they were with Alberta) as they have recorded tree DBH,
# height of pitch tubes, and number of surrounded red attacked trees in the cluster. For 2017-2022 there is no DBH and
# no pitch tube height. There is a vague guesstimate about the number of trees in the cluater, but it's often expressed as ">100",
# meaning so many they couldn't easily be counted.
#
# This paper is the first in a series of two. We mention a "companion paper", which regards the rest of Alberta, a managed
# landscape, which is under provincial jurisdiction, not federal, and includes the commercial pine forest, and so was treated to an
# intense program of removing and burning infested trees. A half billion dollars was spent doing this over the study period 2006-2022.
# The analysis here is structured as to mirror and link to that broader analysis. We hypothesize that whereas (3) predicts (2) predicts (1)
# in the natural setting of the national parks, where no control was undertaken, the same strength of association is not seen
# in the rest of Alberta. (It's there, just weakened.) Specifically: in the rest of Alberta there is a decoupling of r from R
# in the period 2008-2015, and this is a direct result of control effort. The result is an intense outbreak in Jasper
# that did not materialize in the rest of Alberta.
#
# There was an outbreak in Banff that emerged in that same window 2008-2015 as Jasper, but it had a fraction the intensity of Jasper.
# We suspect this is due to the broad valleys in Jasper that are rich in pine, whereas in Banff the valleys are steeper and the pine is higher
# in elevation, and as de la Giroday et al. (2011) reported for British Columbia, low-elevation pine is a key for connecting populations
# in mountainous terrain to get them to erupt and spread. We don't have pine data to go with the elevations, so we will not
# explicitly test this hypothesis.

# additional packages -------------------------------------------------------------------------

## Access databases
library(DBI)
library(odbc)

## statistical analyses
library(mgcv)
library(plotly)
library(reshape2)
library(segmented)

# map locations from Carroll et al. 2017 ------------------------------------------------------

abr_df <- file.path(dataPath, "FRI", "rvaluesQvalues.csv") |>
  read.csv(header = TRUE, na.strings = ".")
abr_sf.latlon <- st_as_sf(abr_df, coords = c("PLOT_LONG", "PLOT_LAT"))
st_crs(abr_sf.latlon) <- latlon

abr_sf <- st_transform(abr_sf.latlon, targetCRS)

## Check whether the r-values plots in the 2017 FRI report are actually outside JNP and BNP.
## We will use this approach later to plot r-value locations in Jasper.
gg_abmpb <- ggplot() +
  geom_sf(data = ab_sf) +
  geom_sf(data = abr_sf, size = 0.5) +
  geom_sf(data = np_banff, col = "blue") +
  geom_sf(data = np_jasper, col = "darkgreen") +
  theme_bw(base_size = 20) +
  annotation_north_arrow(
    location = "bl",
    which_north = "true",
    pad_x = unit(0.25, "in"),
    pad_y = unit(0.25, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  xlab("Longitude") +
  ylab("Latitude")

ggsave(
  file.path(figPath, "carroll_et_al_2017_map_banff_jasper.png"),
  gg_abmpb,
  height = 10,
  width = 7
)

# Part 1. infestation counts/areas for Jasper & Banff ----------------------------------------------
## from Unger, Roke, Thandi & Brett 1999-2022

## Caption:
## Figure 1. Infestation dynamics of mountain pine beetle in Banff and Jasper National Parks, 1999–2021.
## The left y-axis shows the number of infested trees (log scale), and the right y-axis shows the area
## infested in hectares (log scale). Square markers represent counts; circular markers represent
## area estimates. blue is Banff. Pink is Jasper. The vertical dashed line at 2012 marks the
## transition from tree count data to area-based estimates. Note the steep rise in Jasper infestation
## post-2013, contrasting with the more subdued outbreak in Banff.

ABMtnParksMPB <- file.path(dataPath, "Brett", "UngerRokeBrettBanffJasperCountsAreas.txt") |>
  read.table(header = TRUE)

#ABMtnParksMPB <- ABMtnParksMPB |>
#  mutate(
#    ## Adjusted from 123 to 1123 due to likely underestimation at outbreak start
#    Jasperha = ifelse(year == 2013, 1123, Jasperha)
#  )

if (.Platform$OS == "windows") {
  win.graph(height = 5, width = 8)
  par(mar = c(5, 5, 2, 6))
  ## 100 ha is about 2000 mature trees in AB Mtn Parks (20 trees/ha)
  ## Area after 2013 needs to be scaled between 10 ha and 1 000 000 ha
  ## Count before 2013 needs to be scaled between 100 trees to 100 000 trees,
  ## but 100,000 trees is only 50,000 ha, so uncounted tree count after 2013 could be as high as 2,000,000

  ## Banff will be solid black; Jasper white
  ## counts will be squares; areas will be circles

  ## 100 trees to 10 000 000 trees; low end scaled so that white circle perfectly overlaps white square:
  plot(
    ABMtnParksMPB$year,
    log10(ABMtnParksMPB$BanffCount),
    type = "l",
    xlab = "year",
    ylab = "trees infested (log10 count)",
    ylim = c(1, 11.5)
  )
  points(ABMtnParksMPB$year, log10(ABMtnParksMPB$BanffCount), pch = 15) ## black squares
  lines(ABMtnParksMPB$year, log10(ABMtnParksMPB$JasperCount))
  points(ABMtnParksMPB$year, log10(ABMtnParksMPB$JasperCount), pch = 22, bg = "white") ## white squares

  par(new = TRUE)
  plot(
    ABMtnParksMPB$year,
    log10(ABMtnParksMPB$Jasperha),
    axes = FALSE,
    type = "l",
    xlab = "",
    ylab = "",
    ylim = c(1, 6)
  ) ## 10 ha to 1 000 000 ha
  points(ABMtnParksMPB$year, log10(ABMtnParksMPB$Jasperha), pch = 21, bg = "white") ## white circles
  lines(ABMtnParksMPB$year, log10(ABMtnParksMPB$Banffha))
  points(ABMtnParksMPB$year, log10(ABMtnParksMPB$Banffha), pch = 19) ## black circles
  axis(side = 4)
  mtext("area infested (log10 ha)", side = 4, line = 3)

  abline(v = 2012.5, lty = 3)
  text(2005.3, 6, "count infested")
  text(2017.5, 6, "area infested")
  legend(2003.5, 5.7, pch = c(22, 15), c("Jasper", "Banff"))
  legend(2015.5, 1.9, pch = c(21, 19), c("Jasper", "Banff"))
}

## ggplot version
ABMtnParksMPB_long <- ABMtnParksMPB |>
  tidyr::pivot_longer(
    cols = c(BanffCount, JasperCount, Banffha, Jasperha),
    names_to = c("Park", ".value"),
    names_pattern = "(Banff|Jasper)(Count|ha)"
  ) |>
  mutate(Park = as.factor(Park)) |>
  rename(Year = year, Area_ha = ha)

scaleFact <- 23 ## scaling the second y axis

ABMtnParksMPB_plot <- ggplot(ABMtnParksMPB_long, aes(x = Year)) +
  ## Tree count (primary axis)
  geom_line(aes(y = Count, color = Park), linewidth = 1) +
  geom_point(aes(y = Count, fill = Park), shape = 22, size = 3, color = "black", stroke = 0.5) +

  ## Area infested (secondary axis, scaled)
  geom_line(aes(y = Area_ha * scaleFact, color = Park), size = 1) +
  geom_point(
    aes(y = Area_ha * scaleFact, fill = Park),
    shape = 21,
    size = 3,
    color = "black",
    stroke = 0.5
  ) +

  ## Outbreak onset marker
  #geom_vline(xintercept = 2012, linetype = "dashed", linewidth = 1.5, col="gray") +

  ## Log-scaled y-axis with natural tick labels
  scale_y_continuous(
    transform = "log10",
    name = "Trees Infested (count)",
    limits = c(10, 3e7),
    breaks = c(10, 100, 1000, 10000, 1e5, 1e6, 1e7),
    labels = label_number(),

    sec.axis = sec_axis(
      transform = ~ . / scaleFact,
      name = "Area Infested (ha)",
      breaks = c(10, 100, 1000, 10000, 1e5, 1e6, 1e7),
      labels = label_number()
    )
  ) +
  ## Text labels
  geom_text(
    data = data.frame(
      x = c(2004, 2018),
      y = c(2e7, 2e7),
      label = c("     Count Infested", "Area Infested")
    ),
    aes(x = x, y = y, label = label),
    size = 5,
    inherit.aes = FALSE
  ) +
  ## Theme and legend styling
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) +
  scale_x_continuous(limits = c(1998, 2024)) +
  ## Color and fill scales for consistent legend appearance
  scale_fill_manual(values = c("Banff" = "#56B4E9", "Jasper" = "#e75480")) +
  scale_color_manual(values = c("Banff" = "#56B4E9", "Jasper" = "#e75480"))

ggsave(
  file.path(figPath, "AB_mtn_parks_infested_gg.png"),
  ABMtnParksMPB_plot,
  height = 5,
  width = 7
)
ggsave(
  file.path(figPath, "AB_mtn_parks_infested_gg.pdf"),
  ABMtnParksMPB_plot,
  height = 5,
  width = 7
)

## Compute Rt interannual rate of change in area infested Rt=At+1/At for each area
## Note the index Rt is unitless, so we can treat tree counts and areas identically.
compute_rt <- function(x) {
  c(NA, x[-1] / x[-length(x)])
}

## Compute Rt for Banff using area if available, else count
banff_values <- coalesce(ABMtnParksMPB$Banffha, ABMtnParksMPB$BanffCount)
Rt_Banff <- compute_rt(banff_values)

## Compute Rt for Jasper using area if available, else count
jasper_values <- coalesce(ABMtnParksMPB$Jasperha, ABMtnParksMPB$JasperCount)
Rt_Jasper <- compute_rt(jasper_values)

JB.cor<-cor(jasper_values,banff_values,use = "pairwise.complete.obs")
cat("The Jasper-Banff correlation in A/C 1999-2023 is:", JB.cor)

## Combine
JB.Rt <- data.frame(
  year = ABMtnParksMPB$year,
  Rt_Jasper = Rt_Jasper,
  Rt_Banff = Rt_Banff
)

JB.Rt.cor <- cor(JB.Rt$Rt_Jasper, JB.Rt$Rt_Banff, use = "pairwise.complete.obs")
cat("The Jasper-Banff correlation in Rt 1999-2023 is:", JB.Rt.cor)
cor.test(JB.Rt$Rt_Jasper, JB.Rt$Rt_Banff)

# Part 2. r-values for Jasper ----------------------------------------------------------------------

## Processing the source mdb files

mdb_dir <- file.path(dataPath, "Brett", "MtnParksmdb")
mdb_files <- list.files(mdb_dir, pattern = "\\.mdb$", full.names = TRUE)

output_dir <- file.path(mdb_dir, "extracted")
file.path(output_dir, c("source", "site", "tree")) |> fs::dir_create()

if (extract_mdb) {
  stopifnot(.Platform$OS.type == "windows")

  walk(mdb_files, function(mdb) {
    tmp_mdb <- tempfile(fileext = ".mdb")
    file.copy(mdb, tmp_mdb)

    con <- dbConnect(
      odbc::odbc(),
      .connection_string = paste0(
        "Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=",
        tmp_mdb
      )
    )

    db_tbls <- dbListTables(con)

    if ("mpb_trees" %in% db_tbls) {
      write.csv(
        dbReadTable(con, "mpb_trees"),
        file.path(output_dir, "tree", paste0(basename(mdb), "_mpb_trees.csv")),
        row.names = FALSE
      )
    }

    if ("mpb_survey_info" %in% db_tbls) {
      write.csv(
        dbReadTable(con, "mpb_survey_info"),
        file.path(output_dir, "site", paste0(basename(mdb), "_mpb_survey_info.csv")),
        row.names = FALSE
      )
    }

    if ("mpb_site" %in% db_tbls) {
      write.csv(
        dbReadTable(con, "mpb_site"),
        file.path(output_dir, "site", paste0(basename(mdb), "_mpb_site.csv")),
        row.names = FALSE
      )
    }

    file.copy(mdb, file.path(output_dir, "source", basename(mdb)))
    dbDisconnect(con)
    unlink(tmp_mdb)
  })
}

## Joining the site and tree data for 2014-2016 on siteID
tree_dir <- file.path(mdb_dir, "extracted/tree")
site_dir <- file.path(mdb_dir, "extracted/site")

tree_files <- list.files(tree_dir, pattern = "_mpb_trees\\.csv$", full.names = TRUE)
site_files <- list.files(
  site_dir,
  pattern = "_mpb_survey_info\\.csv$|_mpb_site\\.csv$",
  full.names = TRUE
)

## Check for siteID presence
check_siteID <- function(file) {
  df <- readr::read_csv(file, n_max = 100, show_col_types = FALSE)
  data.frame(file = basename(file), has_siteID = "siteID" %in% names(df))
}

tree_check <- bind_rows(lapply(tree_files, check_siteID))
site_check <- bind_rows(lapply(site_files, check_siteID))

print(tree_check)
print(site_check)

## Define columns to keep
site_cols <- c(
  "siteID",
  "beetle_yr",
  "elevation",
  "nbr_infested",
  "plot_lat_dd",
  "plot_long_dd",
  "tot_holes",
  "r_value",
  "survival"
)

tree_cols <- c(
  "siteID",
  "dbh",
  "ht_pitch_tube",
  paste0(
    rep(c("ns1_", "ns2_"), each = 5),
    c("larvae_live", "larvae_dead", "pupae_live", "pupae_dead", "adults_live")[1:5]
  ),
  "ns1_holes",
  "ns2_holes",
  paste0(
    rep(c("ss1_", "ss2_"), each = 5),
    c("larvae_live", "larvae_dead", "pupae_live", "pupae_dead", "adults_live")[1:5]
  ),
  "ss1_holes",
  "ss2_holes"
)

## Function to process one year
process_year <- function(year) {
  base <- paste0("PopForecast_Jasper_BY", year, ".mdb")

  site_file <- file.path(site_dir, paste0(base, "_mpb_site.csv"))
  tree_file <- file.path(tree_dir, paste0(base, "_mpb_trees.csv"))

  site <- readr::read_csv(site_file, show_col_types = FALSE) |> dplyr::select(any_of(site_cols))
  tree <- readr::read_csv(tree_file, show_col_types = FALSE) |> dplyr::select(any_of(tree_cols))

  tree_ids <- readr::read_csv(tree_file, show_col_types = FALSE) |> distinct(siteID)
  site_ids <- readr::read_csv(site_file, show_col_types = FALSE) |> distinct(siteID)

  left_join(tree, site, by = "siteID") |> mutate(beetle_yr = as.integer(year))
}

## Process all years
jasper_rvalues.2014.2016.raw <- bind_rows(lapply(2014:2016, process_year))

## There is a site-tree mismatch in SiteID occurring in 2015 and 2016
jasper_rvalues.2014.2016 <- jasper_rvalues.2014.2016.raw |>
  filter(!is.na(siteID) & !is.na(plot_lat_dd) & !is.na(plot_long_dd))

## check for missing lat/lon
jasper_rvalues.2014.2016 |>
  filter(is.na(plot_lat_dd) | is.na(plot_long_dd))

hist(jasper_rvalues.2014.2016$plot_lat_dd)
hist(jasper_rvalues.2014.2016$plot_long_dd)

## compute our own r-values
jasper_custom_rvalues <- jasper_rvalues.2014.2016 |>
  mutate(
    live_total = rowSums(across(
      c(
        ns1_larvae_live,
        ns2_larvae_live,
        ns1_pupae_live,
        ns2_pupae_live,
        ns1_adults_live,
        ns2_adults_live,
        ss1_larvae_live,
        ss2_larvae_live,
        ss1_pupae_live,
        ss2_pupae_live,
        ss1_adults_live,
        ss2_adults_live
      ),
      ~ replace_na(.x, 0)
    )),

    hole_total = rowSums(across(
      c(
        ns1_holes,
        ns2_holes,
        ss1_holes,
        ss2_holes
      ),
      ~ replace_na(.x, 0)
    )),

    hole_total = rowSums(across(
      c(
        ns1_holes,
        ns2_holes,
        ss1_holes,
        ss2_holes
      ),
      ~ replace_na(.x, 0)
    )) +
      1, ## add 1 to avoid division by zero

    r_tree = live_total / hole_total
  )

## Compare our computed r-values to their reported r-values
ggplot(jasper_custom_rvalues, aes(x = r_value, y = r_tree)) +
  geom_point(alpha = 0.5, color = "darkblue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    x = "Provided r-value (site-level)",
    y = "Custom r-value (tree-level)",
    title = "Comparison of Provided vs Custom r-values"
  ) +
  theme_minimal()

with(jasper_custom_rvalues, cor(r_value, r_tree, use = "complete.obs"))
summary(jasper_custom_rvalues$r_value)
summary(jasper_custom_rvalues$r_tree)

jasper_custom_rvalues |>
  dplyr::select(r_value, r_tree) |>
  pivot_longer(cols = everything(), names_to = "source", values_to = "r") |>
  ggplot(aes(x = r, fill = source)) +
  geom_histogram(position = "dodge", alpha = 0.8, bins = 40) +
  scale_fill_manual(values = c("r_value" = "steelblue", "r_tree" = "darkorange")) +
  labs(
    title = "Distribution of Provided vs Custom r-values",
    x = "r-value",
    y = "Count",
    fill = "Source"
  ) +
  theme_minimal(base_size = 14)

mean(jasper_custom_rvalues$r_tree, na.rm = TRUE)
mean(jasper_custom_rvalues$r_value, na.rm = TRUE)
sd(jasper_custom_rvalues$r_tree, na.rm = TRUE)
sd(jasper_custom_rvalues$r_value, na.rm = TRUE)

## Mean:variance ratios
print("mean-variance ratios:")
print("Custom r-value")
sd(jasper_custom_rvalues$r_tree, na.rm = TRUE) / mean(jasper_custom_rvalues$r_tree, na.rm = TRUE)
print("Provided r-value")
sd(jasper_custom_rvalues$r_value, na.rm = TRUE) / mean(jasper_custom_rvalues$r_value, na.rm = TRUE)

## Build a preliminary gams model of r-values 2014-2016

gam_model.jasper_custom <- gam(
  log(r_tree + 1) ~
    beetle_yr +
      s(plot_long_dd, plot_lat_dd, bs = "gp", k = 16) +
      s(dbh) +
      s(ht_pitch_tube) +
      s(log10(nbr_infested + 1)),
  data = jasper_custom_rvalues,
  method = "REML",
  family = gaussian(link = "identity")
)

summary(gam_model.jasper_custom)
gam.check(gam_model.jasper_custom)
qq.gam(gam_model.jasper_custom, pch = 20)

png(file.path(figPath, "gam_model_jasper_custom.png"), height = 1600, width = 1600)
plot(gam_model.jasper_custom, scheme = 2, pages = 1, all.terms = TRUE)
dev.off()

## processing the source shapefiles 2017-2022 (survey years, beetle years are the year before)
shp_dir <- file.path(dataPath, "Brett", "MtnParksShapefiles")
shp_files <- dir(shp_dir, pattern = "\\.shp$", full.names = TRUE)

## Create output directory for extracted tables
output_dir <- file.path(shp_dir, "extracted") |> fs::dir_create()

walk(shp_files, function(shp) {
  shp_data <- st_read(shp, quiet = TRUE)

  ## Strip geometry and write attribute table to CSV
  attr_table <- st_drop_geometry(shp_data)

  out_file <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(basename(shp)), "_attributes.csv")
  )
  write.csv(attr_table, out_file, row.names = FALSE)
})

shp_csv_files <- list.files(file.path(shp_dir, "extracted"), full.names = TRUE)

## Examine column names, which are not matched
walk(
  shp_csv_files,
  ~ {
    cat("\n---", .x, "---\n")
    print(names(readr::read_csv(.x, n_max = 1, show_col_types = FALSE)))
  }
)

## Survey year to beetle attack year map (file naming conventions switch between mdb and shp)
survey_to_beetle <- c(
  "2018" = 2017,
  "2019" = 2018,
  "2021" = 2020,
  "2022" = 2021
)

## Skip 2017 (already processed)
shp_csv_files <- shp_csv_files[!str_detect(shp_csv_files, "2017")]

## Harmonization function. [Brett refers to the data tech who varied the archiving standard across years 2017-2022.]
read_brett_csv <- function(path) {
  survey_yr <- str_extract(path, "20\\d{2}")
  beetle_yr <- survey_to_beetle[[survey_yr]]

  df <- readr::read_csv(path, show_col_types = FALSE)

  lat_col <- names(df)[str_detect(names(df), regex("^lat$", ignore_case = TRUE))][1]
  lon_col <- names(df)[str_detect(names(df), regex("^lon[g]?$", ignore_case = TRUE))][1]
  r_col <- names(df)[str_detect(names(df), regex("r[_ ]?val|pooled", ignore_case = TRUE))][1]

  df |>
    transmute(
      lat = .data[[lat_col]],
      lon = .data[[lon_col]],
      r_value = .data[[r_col]],
      beetle_yr = beetle_yr
    )
}

## Combine all years
brett_rvalues <- map_dfr(shp_csv_files, read_brett_csv)

## There was one longitude in 2017 that was -177; change it to -117
brett_rvalues <- brett_rvalues |>
  mutate(
    lon = if_else(ceiling(lon) == -177, -117 - lon %% 1, lon)
  )

# dev.new()
hist(brett_rvalues$lat)
# dev.new()
hist(brett_rvalues$lon)
# dev.new()
hist(brett_rvalues$r_value)
# dev.new()
boxplot(brett_rvalues$r_value ~ brett_rvalues$beetle_yr)

asrd_2014.2016_locs <- jasper_custom_rvalues |>
  dplyr::select(lat = plot_lat_dd, lon = plot_long_dd) |>
  mutate(source = "ASRD Tree-Level Data (2014–2016)")

brett_locs <- brett_rvalues |>
  dplyr::select(lat, lon) |>
  mutate(source = "Brett’s Site-Pooled Data (2017–2021)")

Jasper_rvalue_locations <- bind_rows(asrd_2014.2016_locs, brett_locs)

jasper_rvalue_sf <- Jasper_rvalue_locations |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

jasper_rvalue_sf_lambert <- st_transform(jasper_rvalue_sf, crs = 3347)
ab_sf_lambert <- st_transform(ab_sf, crs = 3347)
np_banff_lambert <- st_transform(np_banff, crs = 3347)
np_jasper_lambert <- st_transform(np_jasper, crs = 3347)

rvalues.locs.map <- ggplot() +
  geom_sf(data = ab_sf_lambert) +
  geom_sf(data = np_banff_lambert, col = "blue") +
  geom_sf(data = np_jasper_lambert, col = "darkgreen") +
  geom_sf(data = jasper_rvalue_sf_lambert, aes(color = source), alpha = 0.7, size = 1.5) +
  scale_color_manual(
    values = c(
      "Brett’s Site-Pooled Data (2017–2021)" = "steelblue",
      "ASRD Tree-Level Data (2014–2016)" = "darkorange"
    )
  ) +
  labs(
    title = "Spatial Distribution of Jasper R-value Sampling",
    x = "Longitude",
    y = "Latitude",
    color = "Data Source"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  file.path(figPath, "JNPBNP_r_locs.png"),
  rvalues.locs.map,
  height = 8,
  width = 8
)

## Boxplot of jasper_rvalues.2014.2016 and brett_rvalues

JNPBNP_rvalues <- bind_rows(
  jasper_rvalues.2014.2016 |> dplyr::select(beetle_yr, r_value),
  brett_rvalues |> dplyr::select(beetle_yr, r_value)
)

JNPBNP_rvalues_boxplot <- ggplot(
  JNPBNP_rvalues,
  aes(x = factor(beetle_yr), y = log10(r_value + 1))
) +
  geom_boxplot(fill = "skyblue", color = "black", outlier.color = "red", outlier.size = 2) +
  scale_y_continuous(
    breaks = log10(c(2, 3, 6, 11, 21, 51)), ## r + 1
    labels = c(1, 2, 5, 10, 20, 50), ## original r
    name = "r-value"
  ) +
  labs(
    title = "Distribution of r-values by year (2014–2021)",
    x = "year"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  file.path(figPath, "JNPBNP_rvalues_boxplot.png"),
  JNPBNP_rvalues_boxplot,
  height = 6,
  width = 6
)

## compare r-value in year t with interannual rate of change in area infested Rt=At+1/At

ABMtnParksMPB <- ABMtnParksMPB |>
  mutate(
    combined_area = rowSums(across(c(Banffha, Jasperha)), na.rm = TRUE)
  )
ABMtnParks_area <- ABMtnParksMPB |>
  filter(!is.na(combined_area)) |>
  dplyr::select(year, combined_area)
ABMtnParks_Rtarea <- ABMtnParks_area |>
  mutate(Rt = lead(combined_area) / combined_area)

r_summary <- JNPBNP_rvalues |>
  group_by(beetle_yr) |>
  summarise(mean_r = mean(r_value, na.rm = TRUE))

comparison_df <- left_join(
  r_summary,
  ABMtnParks_Rtarea,
  by = c("beetle_yr" = "year")
)

Rvr.mod <- lm(Rt ~ log10(mean_r + 1), data = comparison_df)
Rvr.mod.sum <- summary(Rvr.mod)
r2 <- format(round(Rvr.mod.sum$r.squared, 2), nsmall = 2)
pval <- formatC(Rvr.mod.sum$coefficients[2, "Pr(>|t|)"], format = "f", digits = 5)
annot_text <- paste0("R² = ", r2, "\np = ", pval)

JNPBNP_Rvsr <- ggplot(comparison_df, aes(x = log10(mean_r + 1), y = Rt)) +
  geom_point(size = 3, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey70", alpha = 0.4) +
  scale_x_continuous(
    breaks = log10(c(1 + 1, 2 + 1, 5 + 1, 10 + 1, 20 + 1)),
    labels = c(1, 2, 5, 10, 20),
    name = expression("Mean r-value")
  ) +
  labs(
    #title = "Relationship Between r-value, r<sub>t</sub>,<br>and Interannual Infestation Rate, R<sub>t</sub>,<br>across Jasper and Banff",
    y = expression(R[t + 1] == (A[t + 1] / A[t]))
  )  +
  geom_text(
    data = comparison_df,
    aes(x = log10(mean_r + 1), y = Rt, label = beetle_yr),
    vjust = -1, size = 4
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_blank(), # optional: removes full box
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

JNPBNP_Rvsr <- JNPBNP_Rvsr +
  annotate("text", x = 0.6, y = 3.2,
           label = annot_text, size = 5, hjust = 0.5, vjust = -0.5)

ggsave(
  file.path(figPath, "JNPBNP_RvsR.png"),
  JNPBNP_Rvsr,
  height = 5,
  width = 5
)

summary(lm(Rt ~ log(mean_r + 1), data = comparison_df))
## marginally significant, but partly because we lumped Jasper and Banff areas
## because they were lumped in the r-values.
## We tried splitting the r-values to match the split areas, but the disaggregation
### brings out an anomaly in the Banff data.

# Part 3. BioSIM simulations of Psurv and CMI for Jasper and Banff -------------------------------------------------

## Convert to sf object
brett_sf <- brett_rvalues |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

## Get elevation (units: meters)
brett_elev <- elevatr::get_elev_point(brett_sf, src = "aws", z = 10)

stopifnot(
  nrow(filter(brett_elev, elevation < 0)) == 0
)

## Combine elevation with original data
brett_rvalues_elev <- brett_elev |>
  st_drop_geometry() |>
  dplyr::select(elevation = elevation) |>
  bind_cols(brett_rvalues)

JNPBNP.locyears <- bind_rows(
  jasper_rvalues.2014.2016 |>
    dplyr::select(lat = plot_lat_dd, lon = plot_long_dd, beetle_yr, elevation),
  brett_rvalues_elev |>
    dplyr::select(lat, lon, beetle_yr, elevation)
)

## input tibble: JNPBNP.locyears
## Columns: lat, lon, beetle_yr, elevation

## Add survey year (BioSIM MPBwk uses beetle year and survey year to define the winter)
JNPBNP.locyears <- JNPBNP.locyears |>
  mutate(survey_yr = beetle_yr + 1)

## Create a unique ID for each location-year combo
JNPBNP.locyears <- JNPBNP.locyears |>
  mutate(loc_id = paste0("loc_", row_number()))

## Reorder and rename columns for clarity
biosim_input <- JNPBNP.locyears |>
  dplyr::select(loc_id, lat, lon, elevation, survey_yr)

## Export to CSV for BioSIM batch processing
## This file can be used in BioSIM's batch mode or uploaded via its web interface
write.csv(biosim_input, file.path(outputPath, "biosim_input_JNPBNP.csv"))

source("R/biosim.R")

f_MPBwkPsurv <- file.path(outputPath, "MPBwkPsurv.csv")
if (!file.exists(f_MPBwkPsurv)) {
  mpb_cold_tol(JNPBNP.locyears) |> write.csv(f_MPBwkPsurv)
}
MPBwkPsurv <- read.csv(f_MPBwkPsurv)

## join MPB winterkill simulation results to main data table, JNPBNP
JNPBNP <- MPBwkPsurv |>
  left_join(mutate(JNPBNP.locyears, row_index = row_number()), by = "row_index")

unique_years <- sort(unique(JNPBNP$survey_yr))

JNPBNP_sf <- st_as_sf(JNPBNP, coords = c("lon", "lat"), crs = 4326)
ab_sf <- st_transform(ab_sf, st_crs(JNPBNP_sf))
st_crs(JNPBNP_sf) == st_crs(ab_sf)

JNPBNP_Psurv_map <- ggplot(JNPBNP_sf |> filter(Year %in% unique_years)) +
  geom_sf(data = ab_sf, fill = NA, color = "black") +
  geom_sf(data = np_banff, color = "blue") +
  geom_sf(data = np_jasper, color = "darkgreen") +
  geom_sf(aes(color = Psurv), size = 1.5) +
  scale_color_viridis_c(option = "C") +
  facet_wrap(~Year, nrow = 2) +
  coord_sf(xlim = c(-125, -110), ylim = c(48, 60), expand = FALSE) +
  theme_minimal() +
  labs(
    title = "Predicted Overwinter Survival (Psurv) by Year",
    subtitle = "Jasper & Banff National Parks",
    color = "Psurv (%)"
  )

ggsave(
  file.path(figPath, "JNPBNP_Psurv_draft.png"),
  JNPBNP_Psurv_map,
  height = 8,
  width = 10
)

## Move from draft map to a more final version

## re-project to Alberta
JNPBNP_proj <- st_transform(JNPBNP_sf, crs = 32611)
ab_proj <- st_transform(ab_sf, crs = 32611)
np_banff_proj <- st_transform(np_banff, crs = 32611)
np_jasper_proj <- st_transform(np_jasper, crs = 32611)

JNPBNP_Psurv_map <- ggplot(JNPBNP_proj |> filter(Year %in% unique_years)) +
  geom_sf(data = ab_proj, fill = NA, color = "black") +
  geom_sf(data = np_banff_proj, color = "blue") +
  geom_sf(data = np_jasper_proj, color = "darkgreen") +
  geom_sf(aes(color = Psurv), size = 1.8) +
  scale_color_viridis_c(
    option = "C",
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    name = "Psurv (%)"
  ) +
  facet_wrap(~Year, nrow = 2) +
  coord_sf(
    xlim = c(300000, 700000),
    ylim = c(5475000, 6000000),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.background = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  labs(
    title = "Predicted Overwinter Survival (Psurv) by Year",
    subtitle = "Jasper & Banff National Parks",
    color = "Psurv (%)"
  )

ggsave(
  file.path(figPath, "JNPBNP_Psurv.png"),
  JNPBNP_Psurv_map,
  height = 8,
  width = 10
)

## boxplot
JNPBNP_Psurv_boxplot <- ggplot(
  JNPBNP_sf |> filter(Year %in% unique_years),
  aes(x = factor(Year), y = Psurv)
) +
  geom_boxplot(fill = "skyblue", color = "darkblue", outlier.color = "red") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  labs(
    title = "Distribution of Predicted Overwinter Survival (Psurv) by Survey Year",
    x = "Year",
    y = "Psurv (%)"
  )

ggsave(
  file.path(figPath, "JNPBNP_Psurv_boxplot.png"),
  JNPBNP_Psurv_boxplot,
  height = 5,
  width = 7
)

## Recall JNPBNP_rvalues was created as follows:
# JNPBNP_rvalues <- bind_rows(
#   jasper_rvalues.2014.2016 |> dplyr::select(beetle_yr, r_value),
#   brett_rvalues |> dplyr::select(beetle_yr, r_value)
# )
## And recall the location-years were created as follows using the same source objects:
# JNPBNP.locyears <- bind_rows(
#   jasper_rvalues.2014.2016 |>
#     dplyr::select(lat = plot_lat_dd, lon = plot_long_dd, beetle_yr, elevation),
#   brett_rvalues_elev |>
#     dplyr::select(lat, lon, beetle_yr, elevation)
#)

## We want to join JNPBNP, which contains Psurv, to JNPBNP_rvalues. which contains nothing but beetle_yr r-values.
## But there is nothing to join them on. So we need to re-create the r-values tibble with some
## simulation point id's so the tables can be joined safely, .i.e not using cbind.

JNPBNP.r <- bind_rows(
  jasper_rvalues.2014.2016 |>
    dplyr::select(beetle_yr, r_value, lat = plot_lat_dd, lon = plot_long_dd, beetle_yr, elevation),
  brett_rvalues_elev |>
    dplyr::select(beetle_yr, r_value, lat = lat, lon = lon, beetle_yr, elevation)
) |>
  mutate(KeyID = paste0("site_", row_number()))

## JNPBNP.r now has the same row structure as JNPBNP
stopifnot(
  nrow(JNPBNP.r) == nrow(JNPBNP),
  all(JNPBNP.r$KeyID %in% JNPBNP$KeyID),
  all.equal(
    as.data.frame(JNPBNP.r |> dplyr::select(lat, lon, beetle_yr) |> distinct()),
    as.data.frame(JNPBNP |> dplyr::select(lat, lon, beetle_yr) |> distinct())
  )
)

## modelling and plotting rt (in JNPBNP.r) on Psurv (in JNPBNP)
JNPBNP.full <- JNPBNP.r |>
  left_join(JNPBNP |> dplyr::select(KeyID, Psurv), by = "KeyID")

#report on r-value means before and after 2017
JNPBNP.full |>
  mutate(period = if_else(beetle_yr %in% 2014:2017, "2014–2017", "Post-2017")) |>
  group_by(period) |>
  summarise(mean_r = mean(r_value, na.rm = TRUE))

#report on r-value means before and after 2017
JNPBNP.full |>
  mutate(period = if_else(beetle_yr %in% 2014:2017, "2014–2017", "Post-2017")) |>
  group_by(period) |>
  summarise(mean_r = mean(Psurv, na.rm = TRUE))

n_unique <- JNPBNP.full |> distinct() |> nrow()
cat("Number of unique location-years:", n_unique)

JNPBNP.mod <- lm(r_value ~ Psurv, data = JNPBNP.full)
JNPBNP.mod.sum<-summary(JNPBNP.mod)

r2<-JNPBNP.mod.sum$r.squared
pval<-pf(JNPBNP.mod.sum$fstatistic["value"],
   JNPBNP.mod.sum$fstatistic["numdf"],
   JNPBNP.mod.sum$fstatistic["dendf"],
   lower.tail = FALSE)
r2 <- format(round(JNPBNP.mod.sum$r.squared, 2), nsmall = 2)
pvalue <- formatC(pval, format = "f", digits = 5)

JNPBNP.rvsPsurv <- ggplot(JNPBNP.full, aes(x = Psurv, y = r_value)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  ) +
  labs(
    #title = "Relationship Between Overwinter Survival (Psurv) and r-value",
    x = "Overwinter Survival (%)",
    y = "r-value"
  )
annot_text <- paste0("R² = ", r2, "\np = ", pvalue)

JNPBNP.rvsPsurv <- JNPBNP.rvsPsurv +
  annotate("text", x = 45, y = 22,
           label = annot_text, size = 5, hjust = 0.5, vjust = -0.5)

ggsave(
  file.path(figPath, "JNPBNP_r_vs_Psurv.png"),
  JNPBNP.rvsPsurv,
  height = 5,
  width = 5
)

## try aggregating the data by year:
JNPBNP.by.year <- JNPBNP.full |>
  group_by(beetle_yr) |>
  summarise(
    mean_r = mean(r_value, na.rm = TRUE),
    mean_Psurv = mean(Psurv, na.rm = TRUE),
    sd_r = sd(r_value, na.rm = TRUE),
    sd_Psurv = sd(Psurv, na.rm = TRUE),

    .groups = "drop"
  )
JNPBNP.year.mod <- lm(mean_r ~ mean_Psurv, data = JNPBNP.by.year)
JNPBNP.year.mod.sum<-summary(JNPBNP.year.mod)

r2 <- format(round(JNPBNP.year.mod.sum$r.squared, 2), nsmall = 2)
pval<-pf(JNPBNP.year.mod.sum$fstatistic["value"],
         JNPBNP.year.mod.sum$fstatistic["numdf"],
         JNPBNP.year.mod.sum$fstatistic["dendf"],
         lower.tail = FALSE)
pvalue <- formatC(pval, format = "f", digits = 5)

JNPBNP.by.year.plot <- ggplot(JNPBNP.by.year, aes(x = mean_Psurv, y = mean_r)) +
  geom_point(size = 3, color = "black") +
  geom_text(aes(label = beetle_yr), vjust = -1, size = 3.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  ) +
  labs(
    #title = "Yearly Aggregated Relationship Between Psurv and r-value",
    x = "Mean Overwinter Survival (%)",
    y = "Mean r-value"
  )

annot_text <- paste0("R² = ", r2, "\np = ", pvalue)

JNPBNP.by.year.plot <- JNPBNP.by.year.plot +
  annotate("text", x = 45, y = 15,
           label = annot_text, size = 5, hjust = 0.5, vjust = -0.5)

ggsave(
  file.path(figPath, "JNPBNP_r_vs_Psurv_yearly.png"),
  JNPBNP.by.year.plot,
  height = 5,
  width = 5
)

## run BioSIM on the Jasper/Banff locations

two_locations <- tibble(
  id = c("Banff_Airport", "Jasper_Airport"),
  lat = c(51.179, 52.990),
  lon = c(-115.570, -118.058),
  elevation = c(1397, 1020)
)

years <- 1998:2023

biosim_input.JNPBNP <- tidyr::crossing(
  two_locations,
  beetle_yr = years
) |>
  arrange(id, beetle_yr)

f_JNPBNP.MPBwkPsurv <- file.path(outputPath, "JNPBNP.MPBwkPsurv.csv")
if (!file.exists(f_JNPBNP.MPBwkPsurv)) {
  mpb_cold_tol(biosim_input.JNPBNP) |> write.csv(f_JNPBNP.MPBwkPsurv)
}
JNPBNP.MPBwkPsurv <- read.csv(f_JNPBNP.MPBwkPsurv)

JNPBNP.MPBwkPsurv <- JNPBNP.MPBwkPsurv |>
  mutate(location = ifelse(Latitude == 51.179, "Banff", "Jasper"))

JNPBNP_1998_2023_Psurv.ts <- ggplot(
  JNPBNP.MPBwkPsurv,
  aes(x = Year - 1, y = Psurv, color = location) #we plot beetle year (the fall year), not the spring year (the year BioSIM reports when output is completed)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Banff" = "#56B4E9", "Jasper" = "#e75480")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  ) +
  labs(
    title = "Simulated Overwinter Survival (Psurv)",
    x = "Beetle Year",
    y = "Overwinter Survival (%)",
    color = "Location"
  ) +
  scale_x_continuous(limits = c(1998, 2024))

ggsave(
  file.path(figPath, "JNPBNP_1998_2023.png"),
  JNPBNP_1998_2023_Psurv.ts,
  height = 4,
  width = 6
)

## compute correlation
twolocs.wide <- JNPBNP.MPBwkPsurv |>
  dplyr::select(Year, location, Psurv) |>
  pivot_wider(names_from = location, values_from = Psurv)

Psurv.cor<-cor(twolocs.wide$Banff, twolocs.wide$Jasper, use = "complete.obs")

cat("The correlation betwen Jasper and Banff in Psurv is:",Psurv.cor)

# Filter for the desired years
twolocs.subset <- twolocs.wide |>
  filter(Year %in% 2012:2017)

# Compute means
mean_Jasper <- mean(twolocs.subset$Jasper, na.rm = TRUE)
mean_Banff  <- mean(twolocs.subset$Banff, na.rm = TRUE)

# Compute difference (Jasper minus Banff)
diff_Psurv <- mean_Jasper - mean_Banff

# Print results
cat("Mean Psurv (Jasper):", mean_Jasper, "\n")
cat("Mean Psurv (Banff):", mean_Banff, "\n")
cat("Difference (Jasper - Banff):", diff_Psurv, "\n")

## CMI
f_JNPBNP.CMI <- file.path(outputPath, "JNPBNPCMI.csv")
if (!file.exists(f_JNPBNP.CMI)) {
  biosim_cmi(biosim_input.JNPBNP) |> write.csv(f_JNPBNP.CMI, row.names = FALSE)
}
JNPBNP.CMI <- read.csv(f_JNPBNP.CMI)

JNPBNP.CMI <- JNPBNP.CMI |>
  mutate(location = ifelse(Latitude == 51.179, "Banff", "Jasper"))

JNPBNP.CMI.plot <- ggplot(JNPBNP.CMI, aes(x = Year, y = CMI, color = location)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Jasper" = "#e75480", "Banff" = "#56B4E9")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  ) +
  labs(
    title = "Climate Moisture Index (CMI)",
    x = "Beetle Year",
    y = "CMI",
    color = "Location"
  ) +
  scale_x_continuous(limits = c(1998, 2024))

ggsave(
  file.path(figPath, "JNPBNP_CMI_ts.png"),
  JNPBNP.CMI.plot,
  height = 4,
  width = 6
)

JNPBNP.CMI |>
  filter(Year %in% 2013:2017) |>
  group_by(location) |>
  summarise(mean_CMI = mean(CMI, na.rm = TRUE)) |>
  summarise(diff_CMI = diff(mean_CMI))

CMI_means <- JNPBNP.CMI |>
  filter(Year %in% 2013:2017) |>
  group_by(location) |>
  summarise(mean_CMI = mean(CMI, na.rm = TRUE))

diff_CMI <- diff(CMI_means$mean_CMI)

jasper.cmi <- subset(JNPBNP.CMI, location == "Jasper")
banff.cmi <- subset(JNPBNP.CMI, location == "Banff")
cor.jb <- cor(banff.cmi$CMI, jasper.cmi$CMI, use = "complete.obs")
cat("The Jasper-Banff CMI correlation is:", cor.jb, "\n")

## Re-visit the r vs Psurv regression relationship (JNPBNP.year.mod) and include CMI (aggregated by year)
## data = JNPBNP.by.year
JNPBNP.by.year
summary(JNPBNP.year.mod)
JNPBNP.by.year.plot

JNPBNP.by.year$beetle_yr

jasper_cmi_subset <- jasper.cmi |>
  dplyr::filter(Year %in% JNPBNP.by.year$beetle_yr) |>
  dplyr::select(Year, CMI)
banff_cmi_subset <- banff.cmi |>
  dplyr::filter(Year %in% JNPBNP.by.year$beetle_yr) |>
  dplyr::select(Year, CMI)

mean_cmi <- tibble(
  Year = jasper_cmi_subset$Year,
  CMI_mean = (jasper_cmi_subset$CMI + banff_cmi_subset$CMI) / 2
)
JNPBNP.by.year <- JNPBNP.by.year |>
  left_join(mean_cmi, by = c("beetle_yr" = "Year"))

JNPBNP.year.mod <- lm(mean_r ~ mean_Psurv + CMI_mean, data = JNPBNP.by.year)
summary(JNPBNP.year.mod)

# Part 4. Modeling Rt directly from Psurv and CMI for Jasper and Banff -------------------------------------------------
## Recall JB.Rt is the rate of change from year t to t+1 in either tree count or area
## JNPBNP.CMI is the CMI object 1999-2024

JB.Rt
head(JNPBNP.CMI)
head(JNPBNP.MPBwkPsurv)

JB.Rt.long <- JB.Rt |>
  pivot_longer(cols = c(Rt_Jasper, Rt_Banff), names_to = "location", values_to = "Rt") |>
  mutate(location = if_else(location == "Rt_Jasper", "Jasper", "Banff"))

## Summarize CMI by year and location
cmi_summary <- JNPBNP.CMI |>
  group_by(Year, location) |>
  summarise(CMI = mean(CMI, na.rm = TRUE), .groups = "drop")

## Summarize Psurv by year and location
psurv_summary <- JNPBNP.MPBwkPsurv |>
  group_by(Year, location) |>
  summarise(Psurv = mean(Psurv, na.rm = TRUE), .groups = "drop")

## Join the two summaries into one tibble
climate_drivers <- left_join(cmi_summary, psurv_summary, by = c("Year", "location")) |>
  rename(year = Year)

Rt_model_data <- JB.Rt.long |>
  left_join(climate_drivers, by = c("year", "location"))

Rt_model_data <- Rt_model_data |>
  arrange(location, year) |>
  group_by(location) |>
  mutate(
    CMI_lag = lag(CMI),
    Psurv_lag = lag(Psurv)
  ) |>
  ungroup()

Rt_lm_pooled <- lm(Rt ~ CMI + CMI_lag + Psurv_lag, data = Rt_model_data)
summary(Rt_lm_pooled)

## try a thresholded model
Rt_model_data_thresh <- Rt_model_data |>
  mutate(CMI_thresh = CMI_lag < -20)

Rt_thresh_model <- lm(Rt ~ CMI_thresh + Psurv_lag, data = Rt_model_data_thresh)
summary(Rt_thresh_model)

Rt_gam <- gam(Rt ~ s(CMI_lag) + Psurv_lag, data = Rt_model_data, method = "REML")
summary(Rt_gam)

Rt_model_data <- Rt_model_data |>
  mutate(location = factor(location))

Rt_gam_interact <- gam(
  Rt ~ s(CMI_lag, by = location) + location + Psurv_lag,
  data = Rt_model_data,
  method = "REML"
)
summary(Rt_gam_interact)

## Come back to the threshold model using segmentation

## Start with a linear model
lm_base <- lm(Rt ~ CMI_lag + Psurv_lag, data = Rt_model_data)

## Fit segmented model with estimated breakpoint in CMI_lag
seg_model <- segmented(lm_base, seg.Z = ~CMI_lag)
summary(seg_model)

## estimates the segmentation at CMI = -12.5

## revisit the threshold model
Rt_model_data_thresh <- Rt_model_data |>
  mutate(CMI_thresh = CMI_lag < -22)

Rt_thresh_model <- lm(Rt ~ CMI_thresh + Psurv_lag, data = Rt_model_data_thresh)
summary(Rt_thresh_model)

grid <- expand.grid(
  CMI_lag = seq(
    min(Rt_model_data$CMI_lag, na.rm = TRUE),
    max(Rt_model_data$CMI_lag, na.rm = TRUE),
    length.out = 100
  ),
  Psurv_lag = seq(
    min(Rt_model_data$Psurv_lag, na.rm = TRUE),
    max(Rt_model_data$Psurv_lag, na.rm = TRUE),
    length.out = 100
  )
)

grid <- grid |>
  mutate(CMI_thresh = CMI_lag < -22)

grid <- grid |>
  mutate(Rt_pred = predict(Rt_thresh_model, newdata = grid))

## Reshape grid into a matrix for z-values
z_matrix <- matrix(
  grid$Rt_pred,
  nrow = length(unique(grid$Psurv_lag)),
  ncol = length(unique(grid$CMI_lag)),
  byrow = TRUE
)

## Extract x and y axes
x_vals <- sort(unique(grid$CMI_lag))
y_vals <- sort(unique(grid$Psurv_lag))

## Create the surface plot
Rt.threhold.model.plot <- plot_ly(x = x_vals, y = y_vals, z = z_matrix, type = "surface") |>
  plotly::layout(
    scene = list(
      xaxis = list(title = "CMI"),
      yaxis = list(title = "Psurv"),
      zaxis = list(title = "Predicted Rt")
    ),
    title = list(text = "Predicted Rt Surface with CMI Threshold at –22")
  )

# Save the plot as a PNG
fig1 <- file.path(figPath, "Rt_threshold_surface.htm")
htmlwidgets::saveWidget(Rt.threhold.model.plot, fig1)

## Use webshot2 (which uses headless Chrome)
fig2 <- file.path(figPath, "Rt_threshold_surface.png")
webshot2::webshot(fig1, fig2, vwidth = 1200, vheight = 900)

# Create surface plot
plot_ly() |>
  add_surface(
    x = x_vals,
    y = y_vals,
    z = z_matrix,
    showscale = FALSE,
    opacity = 0.75
  ) |>
  add_markers(
    data = Rt_model_data_thresh |> filter(location == "Jasper"),
    x = ~CMI_lag,
    y = ~Psurv_lag,
    z = ~Rt,
    marker = list(size = 4, color = '#e75480'),
    name = "Jasper"
  ) |>
  add_markers(
    data = Rt_model_data_thresh |> filter(location == "Banff"),
    x = ~CMI_lag,
    y = ~Psurv_lag,
    z = ~Rt,
    marker = list(size = 4, color = '#56B4E9'),
    name = "Banff"
  ) |>
  layout(
    title = list(text = "Predicted Rt Surface with Observed Data Points"),
    scene = list(
      xaxis = list(title = "CMI"),
      yaxis = list(title = "Psurv"),
      zaxis = list(title = "Rt")
    )
  )

## Figure 1: map of infested areas over DEM

mpb_jb <- st_read(MPB_Jasper_Banff_2012_2023_shp) |>
  st_zm() |>
  st_transform(targetCRS)
glimpse(mpb_jb)

if (plot_all) {
  ggplot() +
    geom_sf(data = st_transform(ab_sf, targetCRS)) + ## Alberta boundary or base layer
    geom_sf(data = np_banff, col = "blue") + ## Banff National Park
    geom_sf(data = np_jasper, col = "darkgreen") + ## Jasper National Park
    geom_sf(data = mpb_jb, fill = "red", alpha = 0.6) + ## MPB polygons in red
    theme_minimal() +
    labs(title = "Mountain Pine Beetle 2012–2023") +
    coord_sf()
}

## Get elevation raster
elev <- elevatr::get_elev_raster(locations = bbox_parks, z = 9, clip = "bbox")

## Convert elevation raster to data frame for ggplot
elev_df <- as.data.frame(raster::rasterToPoints(elev))
colnames(elev_df) <- c("x", "y", "elevation")

mpb_jb$MPB <- "MPB" #add an item for the legend

mpb.map <- ggplot() +
  geom_raster(data = elev_df, aes(x = x, y = y, fill = elevation)) +
  scale_fill_gradientn(colors = terrain.colors(10), name = "Elevation (m)") +
  geom_sf(data = parks, fill = NA, color = "black") +
  geom_sf(data = mpb_jb, aes(color = MPB), fill = "red") +
  scale_color_manual(
    name = NULL,
    values = c("MPB" = "red"),
    guide = guide_legend(override.aes = list(fill = "red"))
  ) +
  theme_minimal() +
  labs(title = "MPB in Jasper & Banff National Parks, 2013–2023") +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_sf()

## add roads and bodies of water to the plot
mpb.map <- mpb.map +
  geom_sf(data = ab_roads_clipped, color = "gray40", size = 0.3) +
  geom_sf(data = ab_hydro_clipped, fill = "lightblue", color = NA) ## TODO: use darker blue?

## Create park label points
townsites_sf <- data.frame(
  name = c("Banff", "Jasper"),
  x = c(-115.57, -118.08),
  y = c(51.176, 52.873)
) |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  st_transform(st_crs(ab_roads))

## add to map (white dot with black outline)
mpb.map <- mpb.map +
  geom_sf(data = townsites_sf, shape = 21, fill = "white", color = "black", size = 3, stroke = 1) +
  geom_sf_text(
    data = townsites_sf,
    aes(label = name),
    nudge_y = 10000, ## adjust as needed for spacing
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  annotation_north_arrow(
    location = "tr", # top right corner
    which_north = "true", # geographic north
    style = north_arrow_fancy_orienteering,
    height = unit(1.5, "cm"),
    width = unit(1.5, "cm")
  ) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  )

## Note: Areas outside the Mountain Parks are shaded to emphasize the scope of the outbreak
## within park boundaries. MPB spread into BC occurred prior to 2013 and is not the focus of this analysis.

## Create a mask polygon
elev_extent <- st_as_sfc(st_bbox(elev))
mask <- st_difference(elev_extent, st_union(parks))

## Add to map
mpb.map <- mpb.map +
  geom_sf(data = mask, fill = "white", alpha = 0.6, color = NA)

ggsave(
  file.path(figPath, "MPB_map_banff_jasper_2013-23.png"),
  mpb.map,
  height = 10,
  width = 10
)

# Generate Final Figures ----------------------------------------------------------------------

## TODO:
library(patchwork)

## Figure 1: map of infested areas over DEM

ggsave(
  file.path(figPath, "Fig1_MPB_map_banff_jasper_2013-23.png"),
  mpb.map,
  height = 12,
  width = 12,
  dpi = 300
)
ggsave(
  file.path(figPath, "Fig1_MPB_map_banff_jasper_2013-23.pdf"),
  mpb.map,
  height = 12,
  width = 12,
  units = "in"
)

## Figure 2: 3-panel time-series
# (a) counts and areas infested 1999-2023
# (b) Psurv 1999-2024
# (c) CMI 1999-2024

#panel(a) ABMtnParksMPB_plot
ABMtnParksMPB_plot.a<- ABMtnParksMPB_plot +
  annotate("text",
           x = min(ABMtnParksMPB_plot$data$Year),
           y = 800000,
           label = "(a)",
           hjust = 0, vjust = 0, size = 6)

#panel(b) JNPBNP_1998_2023_Psurv.ts
#remove the title off panel (b)
JNPBNP_1998_2023_Psurv.ts.b <- JNPBNP_1998_2023_Psurv.ts + labs(title = NULL)
#insert frame label on (b)
JNPBNP_1998_2023_Psurv.ts.b <- JNPBNP_1998_2023_Psurv.ts.b +
  annotate("text", x = min(JNPBNP.MPBwkPsurv$Year) - 1, y = max(JNPBNP.MPBwkPsurv$Psurv),
           label = "(b)", hjust = 0, vjust = 0.5, size = 6)

#panel(c) JNPBNP_CMI
JNPBNP.CMI.plot.c <- JNPBNP.CMI.plot + labs(title = NULL)
JNPBNP.CMI.plot.c<-JNPBNP.CMI.plot.c +
  annotate("text", x = min(JNPBNP.CMI$Year) - 1, y = max(JNPBNP.CMI$CMI),
           label = "(c)", hjust = 0, vjust = 1.2, size = 6)

Fig2_three_panel_plot <-
  ABMtnParksMPB_plot.a /
  (JNPBNP_1998_2023_Psurv.ts.b / JNPBNP.CMI.plot.c + plot_layout(heights = c(3, 3))) +
  plot_layout(heights = c(4, 6))

ggsave(file.path(figPath,"Fig2_JNPBNP_ts.png"), plot = Fig2_three_panel_plot, width = 7, height = 9, units = "in", dpi = 300)
ggsave(file.path(figPath,"Fig2_JNPBNP_ts.pdf"), plot = Fig2_three_panel_plot, width = 7, height = 9, units = "in")

## Figure 3: 2-panel boxplot in time (2014-2022) of
# (a) r-value: JNPBNP_rvalues_boxplot
# (b) Psurv (at plots): JNPBNP_Psurv_boxplot
JNPBNP_rvalues_boxplot
JNPBNP_Psurv_boxplot

custom_theme <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(color = "black"),
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 14),
  plot.title = element_blank(),
)

#panel (a)
JNPBNP_rvalues_boxplot.a <- JNPBNP_rvalues_boxplot + custom_theme +
  xlab("Year of Beetle Attack") +
  annotate("text", x = "2014", y = Inf, label = "(a)", hjust = -0.2, vjust = 2, size = 5)

#panel (b)
JNPBNP_Psurv_boxplot.b <- JNPBNP_Psurv_boxplot +
  custom_theme +
  scale_x_discrete(labels = as.character(c(c(2014:2018),c(2020:2021)))) +
  xlab("Year Winter Starts") +
  ylab("Overwintering Survival (%)") +
  geom_text(data = data.frame(x = 0.75, y = 100),
            aes(x = x, y = y, label = "(b)"),
            inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 5)

Fig3_two_panel_plot <- JNPBNP_rvalues_boxplot.a /
  JNPBNP_Psurv_boxplot.b

ggsave(file.path(figPath,"Fig3_boxplot.png"), plot = Fig3_two_panel_plot, width = 6, height = 6, units = "in", dpi = 300)
ggsave(file.path(figPath,"Fig3_boxplot.pdf"), plot = Fig3_two_panel_plot, width = 6, height = 6, units = "in")

## Figure 4: 3-panel plot of:
# (a) r-value vs Psurv (scatter)              JNPBNP_r_vs_Psurv
# (b) r-value vs. Psurv (aggregated by year)  JNPBNP_r_vs_Psurv_yearly
# (c) Rt vs rt (aggregated by year)           JNPBNP_RvsR

#panel(a)
JNPBNP.rvsPsurv.a<-JNPBNP.rvsPsurv +
  annotate("text", x = 10, y = 30, label = "(a)", hjust = -0.2, vjust = 2, size = 5)

#panel(b)
JNPBNP.by.year.plot.b<-JNPBNP.by.year.plot +
  annotate("text", x = 20, y = 20, label = "(b)", hjust = -0.2, vjust = 2, size = 5)

#panel(c)
JNPBNP_Rvsr.c<-JNPBNP_Rvsr +
  annotate("text", x = 0.15, y = 4, label = "(c)", hjust = -0.2, vjust = 2, size = 5)

Fig4_three_panel_plot <- JNPBNP.rvsPsurv.a /
  JNPBNP.by.year.plot.b /
  JNPBNP_Rvsr.c

ggsave(file.path(figPath,"Fig4_reg.png"), plot = Fig4_three_panel_plot, width = 6, height = 12, units = "in", dpi = 300)
ggsave(file.path(figPath,"Fig4_reg.pdf"), plot = Fig4_three_panel_plot, width = 6, height = 12, units = "in")
