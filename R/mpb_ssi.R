get_SSI <- function(dsn, year) {
  ssi <- sf::st_read(dsn = dsn, layer = paste0("MPB_SSI_", year)) |>
    sf::st_cast("MULTIPOLYGON")

  ssi <- sf::st_make_valid(ssi)
  ## keep only the SSI values and polygon geometries,
  ## and use 'SSI' as the column/field name
  ssi <- ssi |>
    select(matches("^(MPB_SSI|SSI)$"), Shape) |>
    rename(any_of(c(SSI = "MPB_SSI")))

  return(ssi)
}
