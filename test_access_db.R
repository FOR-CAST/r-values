## connecting to Access databases only works on Windows
stopifnot(.Platform$OS.type == "windows")

library(DBI)

connection_string <- paste0(
  "Driver={Microsoft Access Driver (*.mdb, *.accdb)};",
  "DBQ=", normalizePath("./Smoky_Samples taken at DBH.mdb")
)

con <- dbConnect(odbc::odbc(), .connection_string = connection_string, timeout = 10)

## NOTE: you should now be able to  see the connection and browse the structure of the
##  database in the RStudio 'Connections' tab.

dbListTables(con)

dbListFields(con, "mpb_trees")

## read in the table as an R data.frame
mpb_trees <- dbReadTable(con, "mpb_trees")

summary(mpb_trees)


## cleanup
dbDisconnect(con)
