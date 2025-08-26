## TODO: this could be made more efficient by calling BioSIM API in batches,
##       and by unique lat/lon by year, instead of one row at a time

mpb_cold_tol <- function(df) {
  purrr::map_dfr(1:nrow(df), function(i) {
    row <- df[i, ]

    message(
      "Running BioSIM for row ",
      i,
      " of ",
      nrow(df),
      " (lat: ",
      row$lat,
      ", lon: ",
      row$lon,
      ", year: ",
      row$beetle_yr,
      ")"
    )

    tryCatch(
      {
        result <- BioSIM::generateWeather(
          modelNames = "MPB_Cold_Tolerance_Annual",
          fromYr = row$beetle_yr,
          toYr = row$beetle_yr + 1,
          id = paste0("site_", i),
          latDeg = row$lat,
          longDeg = row$lon,
          elevM = row$elevation,
          rep = 1,
          repModel = 1,
          rcp = "RCP45",
          climModel = "RCM4"
        )

        result[["MPB_Cold_Tolerance_Annual"]] |>
          dplyr::mutate(row_index = i)
      },
      error = function(e) {
        message("❌ Row ", i, " failed: ", e$message)
        return(NULL)
      }
    )
  })
}
