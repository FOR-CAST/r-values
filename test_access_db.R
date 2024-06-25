## connecting to Access databases only works on Windows
stopifnot(.Platform$OS.type == "windows")

library(DBI)
library(odbc)

connection_string <- paste0(
  "Driver={Microsoft Access Driver (*.mdb, *.accdb)};",
  "DBQ=", normalizePath("./data/Smoky_Samples taken at DBH.mdb")
)

connection_string <- paste0(
  "Driver={Microsoft Access Driver (*.mdb, *.accdb)};",
  "DBQ=", normalizePath("./data/AB/mdb/SourceData2006to2008/Population Forecast Surveys beetle_yr_2006/Smoky_Samples taken at DBH.mdb")
)

connection_string <- paste0(
  "Driver={Microsoft Access Driver (*.mdb, *.accdb)};",
  "DBQ=", normalizePath("./data/AB/mdb/SourceData2006to2008/Population Forecast Surveys beetle_yr_2008/Smoky r-values 2008/mpb_survey_apr_6_2009.mdb")
)

con <- dbConnect(odbc::odbc(), .connection_string = connection_string, timeout = 10)

## NOTE: you should now be able to  see the connection and browse the structure of the
##  database in the RStudio 'Connections' tab.

dbListTables(con)

dbListFields(con, "mpb_trees")

## read in the table as an R data.frame
mpb_trees <- dbReadTable(con, "mpb_trees")
surv_info <- dbReadTable(con, "mpb_survey_info")
mpb_site <- dbReadTable(con, "mpb_site")

summary(mpb_trees)


## cleanup
dbDisconnect(con)
