# In this script we examine three components of mountain pine beetle population
# dynamics in two Rocky mountain parks, Jasper and Banff, where no control was undertaken.
# There are more data for Jasper than Banff because the outbreak was more intense
# at Banff, so the analysis is asymmetric. The three data components are:
# 1. trends in red tree counts/area (capital "R" is the interannual rate of change in Xt+1/Xt)
# 2. r-values (brood productivity) (lowercase "r")
# 3. predicted winter mortality (1-Psurv) (and drought; CMI) (using BioSIM)
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
# spatial cokriging efficiently, but we can't use that functionality in the web API. So we work around that here, in R,
# by sending a high density grid of points across JNP and BNP to simulate. Maps are smoothed over the lat and lon grid.
# We also run the Climate Moisture Index (CMI) model in BioSIM, as dryness has been fingered as a key determinant of outbreak potential,
# and just as winter temperature is known to fluctuate severely across years, so does drought.
#
# For the beetle years 2014-16 the Jasper r-values data are "rich" (as they were with Alberta) as they have recorded tree DBH,
# height of pitch tubes, and number of surrunded red attacked trees in the cluster. For 2017-2022 there is no DBH and
# no pitch tube height. There is a vague guestimate about the number of trees in the cluater, but it's often expressed as ">100",
# meaning so many they couldn't easily be counted.
#
# This paper is the first in a series of two. We mention a "companion paper", which regards the rest of Alberta, a managed
# landscape, which is under provincial jurisdiction, not federal, and includes the commercial pine forest, and so was treated to an
# intense progam of removing and burning infested trees. A half billion dollars was spent doing this over the study period 2006-2022.
# The analysis here is structured as to mirror and link to that broader analysis. We hypothesize that wherease (3) predicts (2) predicts (1)
# in the natural setting of the national parks, where no control was undertaken, the same strength of association is not seen
# in the rest of Alberta. (It's there, just weakened.) Specifically: in the rest of Alberta there is a decoupling of r from R
# in the period 2008-2015, and this is a direct result of control effort. The result is an intense outbreak in Jasper
# that did not materialize in the rest of Alberta.
#
# There was an outbreak in Banff that emerged in that same window 2008-2015 as Japser, but it had a fraction the intensity of Jasper.
# We suspect this is due to the broad valleys in Jasper that are rich in pine, whereas in Banff the valleys are steeper and the pine is higher
# in elevation, and as de la Giroday et al. (2011) reported for British Columbia, low-elevation pine is a key for connecting populations
# in mountainous terrain to get them to erupt and spread. We don't have pine data to go with the elevations, so we will not
# explicitly test this hypothesis.

# packages ------------------------------------------------------------------------------------

# library(archive)
library(dplyr)
library(ggplot2)
library(ggspatial)
# library(googledrive)
library(sf)
library(terra)
library(scales) #needed for log scale plotting
library(stringr) #needed for wrangling Jasper r-values from shapefiles
library(elevatr)

# setup ---------------------------------------------------------------------------------------

## paths
dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
figPath <- "figures" |> fs::dir_create()
outputPath <- "outputs" |> fs::dir_create()

## set map projection
latlon <- crs("epsg:4326")
targetCRS <- crs(paste(
  "+proj=aea +lat_1=49 +lat_2=67 +lat_0=0 +lon_0=-112 +x_0=0 +y_0=0",
  "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
))


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

# National Parks ------------------------------------------------------------------------------

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

# map locations from Carroll et al. 2017 ------------------------------------------------------

## TODO: use our newly recomputed r-values
abr_df <- read.csv(
  file.path(dataPath, "FRI", "rvaluesQvalues.csv"),
  header = TRUE,
  na.strings = "."
)
abr_sf.latlon <- st_as_sf(abr_df, coords = c("PLOT_LONG", "PLOT_LAT"))
st_crs(abr_sf.latlon) <- latlon

abr_sf <- st_transform(abr_sf.latlon, targetCRS)

#Check whether the r-values plots in the 2017 FRI report are actually outside JNP and BNP.
#We will use this approach later to plot r-value locations in Jasper.
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

# Part 1. infestation counts/areas for Jasper & Banff -------------------------------------------------
## from Unger, Roke, Thandi & Brett 1999-2022

# Caption:
# Figure 1. Infestation dynamics of mountain pine beetle in Banff and Jasper National Parks, 1999–2021.
# The left y-axis shows the number of infested trees (log scale), and the right y-axis shows the area
# infested in hectares (log scale). Square markers represent counts; circular markers represent
# area estimates. blue is Banff. Pink is Jasper. The vertical dashed line at 2012 marks the
# transition from tree count data to area-based estimates. Note the steep rise in Jasper infestation
# post-2013, contrasting with the more subdued outbreak in Banff.

ABMtnParksMPB <- file.path(dataPath, "Brett", "UngerRokeBrettBanffJasperCountsAreas.txt") |>
  read.table(header = TRUE)

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
  points(ABMtnParksMPB$year, log10(ABMtnParksMPB$BanffCount), pch = 15) # black squares
  lines(ABMtnParksMPB$year, log10(ABMtnParksMPB$JasperCount))
  points(ABMtnParksMPB$year, log10(ABMtnParksMPB$JasperCount), pch = 22, bg = "white") # white squares

  par(new = TRUE)
  plot(
    ABMtnParksMPB$year,
    log10(ABMtnParksMPB$Jasperha),
    axes = FALSE,
    type = "l",
    xlab = "",
    ylab = "",
    ylim = c(1, 6)
  ) # 10 ha to 1 000 000 ha
  points(ABMtnParksMPB$year, log10(ABMtnParksMPB$Jasperha), pch = 21, bg = "white") # white circles
  lines(ABMtnParksMPB$year, log10(ABMtnParksMPB$Banffha))
  points(ABMtnParksMPB$year, log10(ABMtnParksMPB$Banffha), pch = 19) # black circles
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

scaleFact <- 16  #scaling the second y axis

ABMtnParksMPB_plot <- ggplot(ABMtnParksMPB_long, aes(x = Year)) +
  # Tree count (primary axis)
  geom_line(aes(y = Count, color = Park), size = 1) +
  geom_point(aes(y = Count, fill = Park), shape = 22, size = 3, color = "black", stroke = 0.5) +

  # Area infested (secondary axis, scaled)
  geom_line(aes(y = Area_ha * scaleFact, color = Park), size = 1) +
  geom_point(aes(y = Area_ha * scaleFact, fill = Park), shape = 21, size = 3, color = "black", stroke = 0.5) +

  # Outbreak onset marker
  geom_vline(xintercept = 2012.5, linetype = "dotted", size = 1.5) +

  # Log-scaled y-axis with natural tick labels
  scale_y_continuous(
    transform = "log10",
    name = "trees infested (count)",
    breaks = c(10, 100, 1000, 10000, 1e5, 1e6, 1e7),
    labels = label_number(),

    sec.axis = sec_axis(
      transform = ~ . / scaleFact,
      name = "area infested (ha)",
      breaks = c(10, 100, 1000, 10000, 1e5, 1e6),
      labels = label_number()
    )
  ) +
  # Text labels
  geom_text(
    data = data.frame(
      x = c(2005.3, 2017.5),
      y = c(1e7, 1e7),
      label = c("count infested", "area infested")
    ),
    aes(x = x, y = y, label = label),
    size = 5,
    inherit.aes = FALSE
  ) +
  # Theme and legend styling
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) +
  # Color and fill scales for consistent legend appearance
  scale_fill_manual(values = c("Banff" = "#56B4E9", "Jasper" = "#e75480")) +
  scale_color_manual(values = c("Banff" = "#56B4E9", "Jasper" = "#e75480"))

ggsave(
  file.path(figPath, "AB_mtn_parks_infested_gg.png"),
  ABMtnParksMPB_plot,
  height = 5,
  width = 7
)

# Part 2. r-values for Jasper -------------------------------------------------

# Processing the source mdb files

library(DBI)
library(odbc)
library(sf)
library(dplyr)
library(fs)
library(purrr)
library(readr)
library(tidyr)

mdb_dir <- file.path(dataPath, "Brett","MtnParksmdb")
mdb_files <- list.files(mdb_dir, pattern = "\\.mdb$", full.names = TRUE)

output_dir <- file.path(mdb_dir, "extracted")
dir_create(output_dir, "source")
dir_create(output_dir, "site")
dir_create(output_dir, "tree")

walk(mdb_files, function(mdb) {
  tmp_mdb <- tempfile(fileext = ".mdb")
  file.copy(mdb, tmp_mdb)

  con <- dbConnect(odbc::odbc(), .connection_string = paste0(
    "Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=", tmp_mdb
  ))

  db_tbls <- dbListTables(con)

  if ("mpb_trees" %in% db_tbls) {
    write.csv(dbReadTable(con, "mpb_trees"),
              file.path(output_dir, "tree", paste0(basename(mdb), "_mpb_trees.csv")),
              row.names = FALSE)
  }

  if ("mpb_survey_info" %in% db_tbls) {
    write.csv(dbReadTable(con, "mpb_survey_info"),
              file.path(output_dir, "site", paste0(basename(mdb), "_mpb_survey_info.csv")),
              row.names = FALSE)
  }

  if ("mpb_site" %in% db_tbls) {
    write.csv(dbReadTable(con, "mpb_site"),
              file.path(output_dir, "site", paste0(basename(mdb), "_mpb_site.csv")),
              row.names = FALSE)
  }

  file.copy(mdb, file.path(output_dir, "source", basename(mdb)))
  dbDisconnect(con)
  unlink(tmp_mdb)
})

#Joining the site and tree data for 2014-2016 on siteID
tree_dir<- file.path(mdb_dir,"extracted/tree")
site_dir<- file.path(mdb_dir,"extracted/site")

tree_files <- list.files(tree_dir, pattern = "_mpb_trees\\.csv$", full.names = TRUE)
site_files <- list.files(site_dir, pattern = "_mpb_survey_info\\.csv$|_mpb_site\\.csv$", full.names = TRUE)

# Check for siteID presence
check_siteID <- function(file) {
  df <- read_csv(file, n_max = 100, show_col_types = FALSE)
  tibble(file = basename(file), has_siteID = "siteID" %in% names(df))
}

tree_check <- bind_rows(lapply(tree_files, check_siteID))
site_check <- bind_rows(lapply(site_files, check_siteID))

print(tree_check)
print(site_check)

# Define columns to keep
site_cols <- c("siteID", "beetle_yr", "elevation", "nbr_infested", "plot_lat_dd",
               "plot_long_dd", "tot_holes", "r_value", "survival")

tree_cols <- c("siteID", "dbh", "ht_pitch_tube",
               paste0(rep(c("ns1_", "ns2_"), each = 5), c("larvae_live", "larvae_dead", "pupae_live", "pupae_dead", "adults_live")[1:5]),
               "ns1_holes", "ns2_holes",
               paste0(rep(c("ss1_", "ss2_"), each = 5), c("larvae_live", "larvae_dead", "pupae_live", "pupae_dead", "adults_live")[1:5]),
               "ss1_holes", "ss2_holes")

# Function to process one year
process_year <- function(year) {
  base <- paste0("PopForecast_Jasper_BY", year, ".mdb")

  site_file <- file.path(site_dir, paste0(base, "_mpb_site.csv"))
  tree_file <- file.path(tree_dir, paste0(base, "_mpb_trees.csv"))

  site <- read_csv(site_file, show_col_types = FALSE) |> dplyr::select(any_of(site_cols))
  tree <- read_csv(tree_file, show_col_types = FALSE) |> dplyr::select(any_of(tree_cols))

  tree_ids <- read_csv(tree_file, show_col_types = FALSE) |> distinct(siteID)
  site_ids <- read_csv(site_file, show_col_types = FALSE) |> distinct(siteID)

  left_join(tree, site, by = "siteID") |> mutate(beetle_yr = year)
}

# Process all years
jasper_rvalues.2014.2016.raw <- bind_rows(lapply(2014:2016, process_year))

#There is a site-tree mismatch in SiteID occurring in 2015 and 2016
jasper_rvalues.2014.2016 <- jasper_rvalues.2014.2016.raw |>
  filter(!is.na(siteID) & !is.na(plot_lat_dd) & !is.na(plot_long_dd))

#check for missing lat/lon
jasper_rvalues.2014.2016 |>
  filter(is.na(plot_lat_dd) | is.na(plot_long_dd))

hist(jasper_rvalues.2014.2016$plot_lat_dd)
hist(jasper_rvalues.2014.2016$plot_long_dd)

#compute our own r-values
jasper_custom_rvalues <- jasper_rvalues.2014.2016 |>
  mutate(
    live_total = rowSums(across(c(
      ns1_larvae_live, ns2_larvae_live,
      ns1_pupae_live, ns2_pupae_live,
      ns1_adults_live, ns2_adults_live,
      ss1_larvae_live, ss2_larvae_live,
      ss1_pupae_live, ss2_pupae_live,
      ss1_adults_live, ss2_adults_live
    ), ~ replace_na(.x, 0))),

    hole_total = rowSums(across(c(
      ns1_holes, ns2_holes, ss1_holes, ss2_holes
    ), ~ replace_na(.x, 0))),

    hole_total = rowSums(across(c(
      ns1_holes, ns2_holes, ss1_holes, ss2_holes
    ), ~ replace_na(.x, 0))) + 1,  # add 1 to avoid division by zero

    r_tree = live_total / hole_total
  )

#Compare our computed r-values to their reported r-values
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

#Mean:variance ratios
print("mean-variance ratios:")
print("Custom r-value")
sd(jasper_custom_rvalues$r_tree, na.rm = TRUE)/mean(jasper_custom_rvalues$r_tree, na.rm = TRUE)
print("Provided r-value")
sd(jasper_custom_rvalues$r_value, na.rm = TRUE)/mean(jasper_custom_rvalues$r_value, na.rm = TRUE)

#Build a preliminary gams model of r-values 2014-2016
library(mgcv)

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

#processing the source shapefiles 2017-2022 (survey years, beetle years are the year before)
shp_dir <- file.path(dataPath, "Brett","MtnParksShapefiles")
shp_files <- dir(shp_dir, pattern = "\\.shp$", full.names = TRUE)

# Create output directory for extracted tables
output_dir <- file.path(shp_dir, "extracted")
dir_create(output_dir)

walk(shp_files, function(shp) {
  shp_data <- st_read(shp, quiet = TRUE)

  # Strip geometry and write attribute table to CSV
  attr_table <- st_drop_geometry(shp_data)

  out_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(shp)), "_attributes.csv"))
  write.csv(attr_table, out_file, row.names = FALSE)
})

shp_csv_files <- list.files(file.path(shp_dir, "extracted"), full.names = TRUE)

#Examine column names, which are not matched
map(shp_csv_files, ~ {
  cat("\n---", .x, "---\n")
  print(names(read_csv(.x, n_max = 1, show_col_types = FALSE)))
})

# Survey year to beetle attack year map (file nameing conventions switch between mdb and shp)
survey_to_beetle <- c(
  "2018" = 2017,
  "2019" = 2018,
  "2021" = 2020,
  "2022" = 2021
)

# Skip 2017 (already processed)
shp_csv_files <- shp_csv_files[!str_detect(shp_csv_files, "2017")]

# Harmonization function. [Brett refers to the data tech who varied the archiving standard across years 2017-2022.]
read_brett_csv <- function(path) {
  survey_yr <- str_extract(path, "20\\d{2}")
  beetle_yr <- survey_to_beetle[[survey_yr]]

  df <- read_csv(path, show_col_types = FALSE)

  lat_col <- names(df)[str_detect(names(df), regex("^lat$", ignore_case = TRUE))][1]
  lon_col <- names(df)[str_detect(names(df), regex("^lon[g]?$", ignore_case = TRUE))][1]
  r_col   <- names(df)[str_detect(names(df), regex("r[_ ]?val|pooled", ignore_case = TRUE))][1]

  df |>
    transmute(
      lat = .data[[lat_col]],
      lon = .data[[lon_col]],
      r_value = .data[[r_col]],
      beetle_yr = beetle_yr
    )
}

# Combine all years
brett_rvalues <- map_dfr(shp_csv_files, read_brett_csv)
dev.new()
hist(brett_rvalues$r_value)
dev.new()
boxplot(brett_rvalues$r_value~brett_rvalues$beetle_yr)

asrd_2014.2016_locs <- jasper_custom_rvalues |>
  dplyr::select(lat = plot_lat_dd, lon = plot_long_dd) |>
  mutate(source = "ASRD Tree-Level Data (2014–2016)")

brett_locs <- brett_rvalues |>
  dplyr::select(lat, lon) |>
  mutate(source = "Brett’s Site-Pooled Data (2017–2021)")

Jasper_rvalue_locations <- bind_rows(asrd_2014.2016_locs, brett_locs)

Jasper_rvalue_locations |>
  filter(lon < -130) |>
  print(width = Inf)

#There was one longitude in 2017 that was -177. I changed it to -117.

library(sf)

jasper_rvalue_sf <- Jasper_rvalue_locations |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

jasper_rvalue_sf_lambert <- st_transform(jasper_rvalue_sf, crs = 3347)
ab_sf_lambert     <- st_transform(ab_sf, crs = 3347)
np_banff_lambert  <- st_transform(np_banff, crs = 3347)
np_jasper_lambert <- st_transform(np_jasper, crs = 3347)

rvalues.locs.map<-ggplot() +
  geom_sf(data = ab_sf_lambert) +
  geom_sf(data = np_banff_lambert, col = "blue") +
  geom_sf(data = np_jasper_lambert, col = "darkgreen") +
  geom_sf(data = jasper_rvalue_sf_lambert, aes(color = source), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c(
    "Brett’s Site-Pooled Data (2017–2021)" = "steelblue",
    "ASRD Tree-Level Data (2014–2016)" = "darkorange"
  )) +
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

#Boxplot of jasper_rvalues.2014.2016 and brett_rvalues

JNPBNP_rvalues <- bind_rows(
  jasper_rvalues.2014.2016 |> dplyr::select(beetle_yr, r_value),
  brett_rvalues |> dplyr::select(beetle_yr, r_value)
)

JNPBNP_rvalues_boxplot<-ggplot(JNPBNP_rvalues, aes(x = factor(beetle_yr), y = log(r_value + 1))) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.color = "red", outlier.size = 2) +
  labs(
    title = "Distribution of log(R + 1) by Year (2014–2021)",
    x = "Year",
    y = "log(R + 1)"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  file.path(figPath, "JNPBNP_rvalues_boxplot.png"),
  JNPBNP_rvalues_boxplot,
  height = 6,
  width = 6
)

# Convert to sf object
brett_sf <- brett_rvalues |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Get elevation (units: meters)
brett_elev <- get_elev_point(brett_sf, src = "aws", z = 10)

# Combine elevation with original data
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


