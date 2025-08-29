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

# infestation counts/areas for Jasper & Banff -------------------------------------------------

## from Unger, Roke, Thandi & Brett 1999-2022

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

scaleFact <- 2

gg <- ggplot(ABMtnParksMPB_long, aes(x = Year, col = Park, fill = Park)) +
  geom_point(aes(y = log10(Count)), shape = 22, size = 3) +
  geom_line(aes(y = log10(Count)), size = 1) +
  geom_point(aes(y = log10(Area_ha) * scaleFact), shape = 21, size = 3) +
  geom_line(aes(y = log10(Area_ha) * scaleFact), size = 1) +
  geom_vline(xintercept = 2012.5, linetype = "dotted", size = 1.5) +
  expand_limits(y = c(0, 11.5)) +
  scale_y_continuous(
    name = "trees infested (log10 count)",
    sec.axis = sec_axis(~ . / scaleFact, name = "area infested (log10 ha)")
  ) +
  geom_text(
    data = data.frame(
      x = c(2005.3, 2017.5),
      y = c(11.5, 11.5),
      label = c("count infested", "area infested")
    ),
    aes(x = x, y = y, label = label),
    size = 5,
    inherit.aes = FALSE
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

## TODO: count value for Jasper @ 2013 ends up separated from other count data on left side of fig.

ggsave(file.path(figPath, "AB_mtn_parks_infested_gg.png"), gg)
