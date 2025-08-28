# additional packages -------------------------------------------------------------------------

library(car)
library(mgcv)

# data and model exploration ------------------------------------------------------------------

## split dataset based on pivot year
pivot_year <- 2015

abr.early <- all_data_df |> filter(beetle_yr <= pivot_year)
abr.late <- all_data_df |> filter(beetle_yr > pivot_year)

sapply(
  abr.early[, c(
    "r",
    "beetle_yr",
    "plot_lat_dd_copy",
    "plot_lon_dd_copy",
    "nbr_infested",
    "dbh",
    "ht_pitch_tube"
  )],
  function(x) sum(is.na(x))
)
sapply(
  abr.late[, c(
    "r",
    "beetle_yr",
    "plot_lat_dd_copy",
    "plot_lon_dd_copy",
    "nbr_infested",
    "dbh",
    "ht_pitch_tube"
  )],
  function(x) sum(is.na(x))
)

### try linear regression ---------------------------------------------------------------------

if (FALSE) {
  abr.lm.early <- lm(
    log10(r + 1) ~
      beetle_yr +
        plot_lat_dd_copy +
        plot_lon_dd_copy +
        Q +
        log10(nbr_infested + 1) * dbh * ht_pitch_tube * SSI_2008,
    data = abr.early,
    na.action = na.omit
  )
  summary(abr.lm.early)

  abr.lm.late <- lm(
    log10(r + 1) ~
      beetle_yr +
        plot_lat_dd_copy +
        plot_lon_dd_copy +
        log10(nbr_infested + 1) * dbh * ht_pitch_tube * SSI_2023,
    data = abr.late,
    na.action = na.omit
  )
  summary(abr.lm.late)

  abr.lm.late <- lm(
    log10(r + 1) ~
      beetle_yr +
        plot_lat_dd_copy +
        plot_lon_dd_copy +
        Q +
        log10(nbr_infested + 1) * dbh * ht_pitch_tube * SSI_2023,
    data = abr.late,
    na.action = na.omit
  )
  summary(abr.lm.late)

  ## Examine variance inflation factor associated with multi-colinearity
  car::vif(abr.lm.early)
  car::vif(abr.lm.late)

  ## We are drowning in multi-colinearity due to SSI and Q and lat/lon.
  ## Try dropping lat/lon

  abr.lm.early <- lm(
    log10(r + 1) ~ beetle_yr + Q + log10(nbr_infested + 1) * dbh * ht_pitch_tube * SSI_2008,
    data = abr.early,
    na.action = na.omit
  )
  summary(abr.lm.early)

  abr.lm.late <- lm(
    log10(r + 1) ~ beetle_yr + Q + log10(nbr_infested + 1) * dbh * ht_pitch_tube * SSI_2016,
    data = abr.late,
    na.action = na.omit
  )
  summary(abr.lm.late)

  ## We are still drowning in multi-colinearity due to SSI and Q.
  ## Try dropping Q.
  abr.lm.early <- lm(
    log10(r + 1) ~
      beetle_yr +
        log10(nbr_infested + 1) * dbh * ht_pitch_tube * SSI_2008,
    data = abr.early,
    na.action = na.omit
  )
  summary(abr.lm.early)

  abr.lm.late <- lm(
    log10(r + 1) ~
      beetle_yr +
        log10(nbr_infested + 1) * dbh * ht_pitch_tube * SSI_2023,
    data = abr.late,
    na.action = na.omit
  )
  summary(abr.lm.late)

  ## We are still drowning in multi-colinearity due to 4-way interaction.
  ## Try dropping SSI.
  abr.lm.early <- lm(
    log10(r + 1) ~
      beetle_yr +
        log10(nbr_infested + 1) * dbh * ht_pitch_tube,
    data = abr.early,
    na.action = na.omit
  )
  summary(abr.lm.early)

  abr.lm.late <- lm(
    log10(r + 1) ~
      beetle_yr +
        log10(nbr_infested + 1) * dbh * ht_pitch_tube,
    data = abr.late,
    na.action = na.omit
  )
  summary(abr.lm.late)

  ## What if SSI affects immigration behaviour, not beetle pressure/colonization success?
  ## Try including SSI not as interaction.
  abr.lm.early <- lm(
    log10(r + 1) ~ beetle_yr + SSI_2008 + log10(nbr_infested + 1) * dbh * ht_pitch_tube,
    data = abr.early,
    na.action = na.omit
  )
  summary(abr.lm.early)

  abr.lm.late <- lm(
    log10(r + 1) ~ beetle_yr + SSI_2023 + log10(nbr_infested + 1) * dbh * ht_pitch_tube,
    data = abr.late,
    na.action = na.omit
  )
  summary(abr.lm.late)

  ## SSI ns as main effect. So remove completely and try Q.
  abr.lm.early <- lm(
    log10(r + 1) ~ beetle_yr + Q + log10(nbr_infested + 1) * dbh * ht_pitch_tube,
    data = abr.early,
    na.action = na.omit
  )
  summary(abr.lm.early)

  abr.lm.late <- lm(
    log10(r + 1) ~ beetle_yr + Q + log10(nbr_infested + 1) * dbh * ht_pitch_tube,
    data = abr.late,
    na.action = na.omit
  )
  summary(abr.lm.late)

  car::vif(abr.lm.early)
  car::vif(abr.lm.late)
}

## Even with the simple model the VIFs are very high, especially on HT and DBH;
## A shame because the DBH and HT distributions are nice.
if (plot_all) {
  dev.new()
  HT.hist <- hist(all_data_df$ht_pitch_tube)
  dev.new()
  DBH.hist <- hist(all_data_df$dbh)
}

### try non-linear models ---------------------------------------------------------------------

if (FALSE) {
  gam_model.e <- gam(
    r ~
      s(dbh) +
        s(ht_pitch_tube) +
        s(log10(nbr_infested + 1)) +
        s(Q, bs = "gp") +
        s(SSI_2008, bs = "gp") +
        s(plot_lon_dd_copy, plot_lat_dd_copy, bs = "gp") +
        beetle_yr,
    data = abr.early
  )
  summary(gam_model.e)
  plot(gam_model.e, scheme = 2, pages = 1, all.terms = TRUE)

  gam_model.l <- gam(
    r ~
      s(dbh) +
        s(ht_pitch_tube) +
        s(log10(nbr_infested + 1)) +
        s(Q, bs = "gp") +
        s(SSI_2023, bs = "gp") +
        s(plot_lon_dd_copy, plot_lat_dd_copy, bs = "gp") +
        beetle_yr,
    data = abr.late
  )
  summary(gam_model.l)
  plot(gam_model.l, scheme = 2, pages = 1, all.terms = TRUE)
}

## with winterkill

if (FALSE) {
  gam_model <- gam(
    Psurv_prop ~ s(Tmin, bs = "gp"),
    data = site_year__MPBwk_results,
    family = binomial(link = "logit")
  )

  gam_model <- gam(
    Psurv ~ s(Tmin, bs = "gp"),
    data = site_year__MPBwk_results
  )
  summary(gam_model)

  gam.check(gam_model)
  plot(gam_model, residuals = TRUE)
}

## End of testing initial models
