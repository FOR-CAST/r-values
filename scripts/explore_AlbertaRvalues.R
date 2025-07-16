##
## superseded by scripts/test_access_db.R
##

## try using RODBC package ----------------------------------------------------

library(RODBC)
library(DBI)

cwd <- "C:\\Users\\bcooke\\Documents\\Barry\\Data\\MPB\\R_value_BY2006to2019\\mdb"

filedir.0608 <- "C:\\Users\\bcooke\\Documents\\Barry\\Data\\MPB\\R_value_BY2006to2019\\mdb\\SourceData2006to2008\\Population Forecast Surveys beetle_yr_2006\\"
filedir.0608 <- "C:/Users/bcooke/Documents/Barry/Data/MPB/mdb/Rvalues0608/"

smokydata <- "Smoky0608.mdb"

smoky0608 <- paste(filedir.0608, smokydata, sep = "")

smoky0608.db <- odbcConnectAccess(smoky0608) ## fails with 64-bit Windows

driver.cmd <- paste("Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=", smoky0608, sep = "")

smoky0608.db <- odbcDriverConnect(driver.cmd)

## try using odbc package -----------------------------------------------------

library(DBI)
library(odbc)

smoky0608.db <- odbcConnectAccess(smoky0608) # fails with 64-bit windows
db_file_path <- filedir.0608

connect_to_access_dbi <- function(db_file_path) {
  require(DBI)
  ## make sure that the file exists before attempting to connect
  if (!file.exists(db_file_path)) {
    stop("DB file does not exist at ", db_file_path)
  }
  ## Assemble connection strings
  dbq_string <- paste0("DBQ=", db_file_path)
  driver_string <- "Driver={Microsoft Access Driver (*.mdb, *.accdb)};"
  db_connect_string <- paste0(driver_string, dbq_string)

  myconn <- dbConnect(odbc::odbc(), .connection_string = db_connect_string)
  return(myconn)
}
smoky0608.db <- connect_to_access_dbi(filedir.0608)

library(Hmisc)
mdb_path <- smoky0608
mdb.get(mdb_path, tables = TRUE)
