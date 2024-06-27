# packages ------------------------------------------------------------------------------------

library(archive)
library(ggplot2)
library(ggspatial)
library(sf)

# setup ---------------------------------------------------------------------------------------

## paths
dataPath <- normalizePath("./data", mustWork = FALSE)
figPath <- "figures"
outputPath <- "outputs"

if (!dir.exists(dataPath)) dir.create(dataPath)
if (!dir.exists(figPath)) dir.create(figPath)
if (!dir.exists(outputPath)) dir.create(outputPath)

## set map projection
latlon <- crs("epsg:4326")
targetCRS <- crs(paste(
  "+proj=aea +lat_1=49 +lat_2=67 +lat_0=0 +lon_0=-112 +x_0=0 +y_0=0",
  "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
))

# National Parks ------------------------------------------------------------------------------

## Banff and Jasper National Parks
## see: https://open.canada.ca/data/dataset/e1f0c975-f40c-4313-9be2-beb951e35f4e
np_url <- "https://opendata.arcgis.com/datasets/0fb235ee5bb34e51a825add061dd1a21_0.zip"
np_zip <- file.path(dataPath, basename(np_url))
np_file <- "vw_Places_Public_lieux_public_APCA"
np_fext <- c("cpg", "dbf", "prj", "shp", "shx")

if (!file.exists(np_zip)) {
  download.file(url = np_url, destfile = np_zip)
}

if (!all(file.exists(file.path(dataPath, paste0(np_file, ".", np_fext))))) {
  archive::archive_extract(np_zip, dataPath)
}

natl_prks.latlon <- st_read(file.path(dataPath, paste0(np_file, ".shp")))
natl_prks <- st_transform(natl_prks.latlon, targetCRS)

np_banff.latlon <- natl_prks.latlon[natl_prks.latlon$DESC_EN == "Banff National Park of Canada", ]
np_banff <- natl_prks[natl_prks$DESC_EN == "Banff National Park of Canada", ]

np_jasper.latlon[natl_prks.latlon$DESC_EN == "Jasper National Park of Canada", ]
np_jasper <- natl_prks[natl_prks$DESC_EN == "Jasper National Park of Canada", ]

# map locations from Carroll et al. 2017 ------------------------------------------------------

abr_df <- read.csv(file.path(dataPath, "FRI", "rvaluesQvalues.csv"), header = TRUE, na.strings = ".")
abr_sf.latlon <- st_as_sf(abr_df, coords = c("PLOT_LONG", "PLOT_LAT"))
st_crs(abr_sf.latlon) <- latlon

abr_sf <- st_transform(abr_sf.latlon, targetCRS)

gg_abmpb <- ggplot() +
  geom_sf(data = ab) +
  geom_sf(data = abr_sf, size = 0.5) +
  geom_sf(data = np_banff, col = "blue") +
  geom_sf(data = np_jasper, col = "darkgreen") +
  theme_bw(base_size = 20) +
  annotation_north_arrow(
    location = "bl", which_north = "true",
    pad_x = unit(0.25, "in"), pad_y = unit(0.25, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  xlab("Longitude") +
  ylab("Latitude")

ggsave(file.path(figPath, "carroll_et_al_2017_map_banff_jasper.png"), gg_abmpb, height = 10, width = 7)
