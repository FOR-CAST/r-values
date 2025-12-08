## Kamloops, BC coordinates: 50.6745° N, 120.3273° W; 345 m elevation

lat <- 50.6745
lon <- -120.3273
elev <- 345

CMI_monthly <- BioSIM::generateWeather(
  modelNames = "Climate_Moisture_Index_Monthly",
  fromYr = 2000,
  toYr = 2005,
  id = "Kamloops",
  latDeg = lat,
  longDeg = lon,
  elevM = elev,
  rep = 1,
  repModel = 1,
  rcp = "RCP45",
  climModel = "RCM4"
)

CMI_annual <- BioSIM::generateWeather(
  modelNames = "Climate_Mosture_Index_Annual", ## typo in BioSIM package!
  fromYr = 2000,
  toYr = 2005,
  id = "Kamloops",
  latDeg = lat,
  longDeg = lon,
  elevM = elev,
  rep = 1,
  repModel = 1,
  rcp = "RCP45",
  climModel = "RCM4"
)
