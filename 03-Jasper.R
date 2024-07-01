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

# infestation counts/areas for Jasper & Banff --------------(from Unger, Roke, Thandi & Brett 2003-2022) --------------------------------------

ABMtnParksMPB <- read.table(file.path(dataPath, "Brett", "UngerRokeBrettBanffJasperCountsAreas.txt"), header = TRUE)

win.graph(height=5,width=8)
par(mar=c(5,5,2,6))
# 100 ha is about 2000 mature trees in AB Mtn Parks (20 trees/ha)
#Area after 2013 needs to be scaled between 10 ha and 1 000 000 ha
#Count before 2013 needs to be scaled between 100 trees to 100 000 trees, but 100 000 trees is only 50 000 ha, so uncounted tree count after 2013 could be as high as 2 000 000

#Banff will be solid black; Jasper white
#counts will be squares; areas will be circles

plot(ABMtnParksMPB$year,log10(ABMtnParksMPB$BanffCount),type="l",xlab="year",ylab="trees infested (log10 count)",ylim=c(1,11.5)) # 100 trees to 10 000 000 trees; low end scaled so that white circle perfectly overlaps white square
points(ABMtnParksMPB$year,log10(ABMtnParksMPB$BanffCount),pch=15) #black squares
lines(ABMtnParksMPB$year,log10(ABMtnParksMPB$JasperCount))
points(ABMtnParksMPB$year,log10(ABMtnParksMPB$JasperCount),pch=22,bg="white") #white squares

par(new=TRUE)
plot(ABMtnParksMPB$year,log10(ABMtnParksMPB$Jasperha),axes=F,type="l",ylim=c(1,6),xlab="",ylab="") # 10 ha to 1 000 000 ha
points(ABMtnParksMPB$year,log10(ABMtnParksMPB$Jasperha),pch=21,bg="white") #white circles
lines(ABMtnParksMPB$year,log10(ABMtnParksMPB$Banffha))
points(ABMtnParksMPB$year,log10(ABMtnParksMPB$Banffha),pch=19) #black circles
axis(side=4)
mtext("area infested (log10 ha)", side = 4, line = 3)

abline(v=2012.5,lty=3)
text(2005.3,6,"count infested")
text(2017.5,6,"area infested")
legend(2003.5,5.7,pch=c(22,15),c("Jasper","Banff"))
legend(2015.5,1.9,pch=c(21,19),c("Jasper","Banff"))

