## based on version from achubaty/mpbPine module
## <https://github.com/achubaty/mpbPine/blob/f9a160099b381e8b53289f1d5b251404d8d53d41/mpbPine.R#L293-L330>
prepInputs_ABPine <- function(url, rasterToMatch, layerNames, ...) {
  reproducible::prepInputs(
    url = url,
    targetFile = "AB_PineVolumes_Lambert.gdb",
    rasterizeMask = rasterizeMask,
    layers = layerNames,
    rasterToMatch = rasterToMatch,
    useCache = FALSE,
    fun = rasterizeMask(targetFile = targetFile, layer = layers, rasterToMatch = rasterToMatch),
    ...
  ) |>
    rastTimes10()
}

## based on version from achubaty/mpbPine module
## <https://github.com/achubaty/mpbPine/blob/f9a160099b381e8b53289f1d5b251404d8d53d41/mpbPine.R#L333-L344>
rasterizeMask <- function(targetFile, rasterToMatch, layer) {
  sf::st_read(targetFile, layer = layer) |>
    terra::rasterize(terra::rast(rasterToMatch), field = "PCT_P") |>
    maskInputs(rasterToMatch = rasterToMatch, studyArea = NULL, maskWithRTM = TRUE)
}

## based on version from achubaty/mpbPine module
## <https://github.com/achubaty/mpbPine/blob/f9a160099b381e8b53289f1d5b251404d8d53d41/mpbPine.R#L340-L344>
rastTimes10 <- function(ras) {
  ras * 10
}
