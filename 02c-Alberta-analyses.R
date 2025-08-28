# additional packages -------------------------------------------------------------------------

library(mgcv)

## build model now with Psurv and/or Tmin -----------------------------------------------------

### split early and late beetle years ---------------------------------------------------------

if (FALSE) {
  abr.early <- all_data_df_join_CMI |> filter(beetle_yr <= pivot_year)
  abr.late <- all_data_df_join_CMI |> filter(beetle_yr > pivot_year)

  gam_model.e <- gam(
    r ~
      beetle_yr +
        s(dbh) +
        s(ht_pitch_tube) +
        s(log10(nbr_infested + 1)) +
        s(Q, bs = "gp") +
        s(asin(sqrt(Q)), bs = "gp") +
        s(lon, lat, bs = "gp") +
        s(Tmin, bs = "gp") +
        s(Psurv, bs = "gp") +
        s(PineVol, bs = "gp"),
    data = abr.early,
    method = "REML",
    family = gaussian(link = "identity")
  )
  summary(gam_model.e)

  gam.check(gam_model.e)

  # dev.new()
  plot(gam_model.e, scheme = 2, pages = 1, all.terms = TRUE)

  gam_model.l <- gam(
    r ~
      beetle_yr +
        s(dbh) +
        s(ht_pitch_tube) +
        s(log10(nbr_infested + 1)) +
        s(asin(sqrt(Q)), bs = "gp") +
        s(SSI_2023, bs = "gp") +
        s(lon, lat, bs = "gp") +
        s(Tmin, bs = "gp") +
        s(Psurv, bs = "gp") +
        s(PineVol, bs = "gp"),
    data = abr.late,
    method = "REML",
    family = gaussian(link = "identity")
  )
  summary(gam_model.l)

  gam.check(gam_model.l)

  # dev.new()
  plot(gam_model.l, scheme = 2, pages = 1, all.terms = TRUE)
}

### use full dataset --------------------------------------------------------------------------

## below we doubled k' based on initial gam.check of model using default k-values;
## then doubled k-values for Psurv and PineVol again.
## (see ?mgcv::gam.check and ?mgcv::choose.k)

gam_model.all <- gam(
  r ~
    s(lon, lat, bs = "gp", k = 64) +
      s(beetle_yr) +
      s(dbh) +
      s(ht_pitch_tube) +
      s(log10(nbr_infested + 1)) +
      s(CMI, bs = "gp", k = 22) +
      s(asin(sqrt(Q)), bs = "gp", k = 22) +
      s(PineVol, bs = "gp", k = 44) +
      s(Psurv, bs = "gp", k = 44) +
      s(SSI_2016, bs = "gp", k = 22) +
      s(Tmin, bs = "gp", k = 22),
  data = all_data_df_join_CMI,
  method = "REML",
  family = gaussian(link = "identity")
)
summary(gam_model.all)

gam.check(gam_model.all)

qq.gam(gam_model.all, pch = 20)

# dev.new()
plot(gam_model.all, scheme = 2, pages = 1, all.terms = TRUE)

png(file.path(figPath, "gam_model_all.png"), height = 1600, width = 1600)
plot(gam_model.all, scheme = 2, pages = 1, all.terms = TRUE)
dev.off()

### -------------------------------------------------------------------------------------------

# ssi_crs <- get_SSI(dsn = ssi_gdb, year = 2023) |> sf::st_crs() ## EPSG:3400
ssi_crs <- sf::st_crs(3400)

## bring the geometry back into the sf
all_data_sf <- all_data_df_join_CMI |>
  filter(!is.na(lon) & !is.na(lat)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_make_valid() |>
  st_transform(st_crs(ssi_crs))

# dev.new()
ggplot(all_data_sf) +
  geom_sf(aes(color = log10(r + 1)), size = 1.1, alpha = 0.7) +
  geom_sf(data = ab_sf, fill = NA) +
  facet_wrap(~beetle_yr, ncol = 7, nrow = 2) +
  scale_color_viridis_c(option = "plasma", name = "log₁₀(r + 1)") +
  theme_minimal() +
  labs(
    title = "Spatial Distribution of r-values by Year",
    subtitle = "Log-scaled to handle extreme values",
    caption = "Each point represents a tree-level estimate"
  )

yearly_summary <- all_data_df_join_CMI |>
  group_by(beetle_yr) |>
  summarise(
    mean_r_log = mean(log10(r + 1), na.rm = TRUE),
    se_r_log = sd(log10(r + 1), na.rm = TRUE) / sqrt(n()),
    mean_Psurv = mean(Psurv, na.rm = TRUE),
    se_Psurv = sd(Psurv, na.rm = TRUE) / sqrt(n()),
    mean_Tmin = mean(Tmin, na.rm = TRUE),
    se_Tmin = sd(Tmin, na.rm = TRUE) / sqrt(n())
  )

plot_df <- yearly_summary |>
  pivot_longer(
    cols = c(mean_r_log, se_r_log, mean_Psurv, se_Psurv, mean_Tmin, se_Tmin),
    names_to = "metric",
    values_to = "value"
  )

stats::cor(yearly_summary$mean_r_log, yearly_summary$mean_Psurv)

stats::cor(yearly_summary$mean_r_log, yearly_summary$mean_Tmin)

# dev.new()
ggplot(yearly_summary, aes(x = beetle_yr)) +
  ## Psurv line and ribbon
  geom_ribbon(
    aes(ymin = mean_Psurv - se_Psurv, ymax = mean_Psurv + se_Psurv, fill = "Psurv (%)"),
    alpha = 0.2
  ) +
  geom_line(aes(y = mean_Psurv, color = "Psurv (%)"), size = 1.2) +
  geom_point(aes(y = mean_Psurv, color = "Psurv (%)"), size = 2) +

  ## log(r) line and ribbon (scaled)
  geom_ribbon(
    aes(
      ymin = (mean_r_log - se_r_log) * 200,
      ymax = (mean_r_log + se_r_log) * 200,
      fill = "log₁₀(r + 1)"
    ),
    alpha = 0.2
  ) +
  geom_line(aes(y = mean_r_log * 200, color = "log₁₀(r + 1)"), size = 1.2) +
  geom_point(aes(y = mean_r_log * 200, color = "log₁₀(r + 1)"), size = 2) +

  ## axes and legend
  scale_y_continuous(
    name = "Psurv (%)",
    limits = c(0, 100),
    sec.axis = sec_axis(~ . / 200, name = "log₁₀(r + 1)")
  ) +
  scale_color_manual(values = c("Psurv (%)" = "#00BFC4", "log₁₀(r + 1)" = "#F8766D")) +
  scale_fill_manual(values = c("Psurv (%)" = "#00BFC4", "log₁₀(r + 1)" = "#F8766D")) +
  theme_minimal() +
  labs(
    title = "Mean Psurv and log-scaled r-values Over Time",
    subtitle = "Dual-axis plot with shaded error ribbons",
    x = "Beetle Year",
    color = "Metric",
    fill = "Metric",
    caption = "Psurv shown as percent; r-values log-transformed and scaled for visibility"
  )

# dev.new()
y_scale <- 50
y_shift <- -50
ggplot(yearly_summary, aes(x = beetle_yr)) +
  ## Tmin line and ribbon
  geom_ribbon(
    aes(ymin = mean_Tmin - se_Tmin, ymax = mean_Tmin + se_Tmin, fill = "Tmin (°C)"),
    alpha = 0.2
  ) +
  geom_line(aes(y = mean_Tmin, color = "Tmin (°C)"), size = 1.2) +
  geom_point(aes(y = mean_Tmin, color = "Tmin (°C)"), size = 2) +

  ## log(r) line and ribbon (scaled)
  geom_ribbon(
    aes(
      ymin = (mean_r_log - se_r_log) * y_scale + y_shift,
      ymax = (mean_r_log + se_r_log) * y_scale + y_shift,
      fill = "log₁₀(r + 1)"
    ),
    alpha = 0.2
  ) +
  geom_line(aes(y = mean_r_log * y_scale + y_shift, color = "log₁₀(r + 1)"), size = 1.2) +
  geom_point(aes(y = mean_r_log * y_scale + y_shift, color = "log₁₀(r + 1)"), size = 2) +

  ## axes and legend
  scale_y_continuous(
    name = "Tmin (°C)",
    limits = c(-50, -20),
    sec.axis = sec_axis(~ (. - y_shift) / y_scale, name = "log₁₀(r + 1)")
  ) +
  scale_color_manual(values = c("Tmin (°C)" = "#00BFC4", "log₁₀(r + 1)" = "#F8766D")) +
  scale_fill_manual(values = c("Tmin (°C)" = "#00BFC4", "log₁₀(r + 1)" = "#F8766D")) +
  theme_minimal() +
  labs(
    title = "Mean Tmin and log-scaled r-values Over Time",
    subtitle = "Dual-axis plot with shaded error ribbons",
    x = "Beetle Year",
    color = "Metric",
    fill = "Metric",
    caption = "Tmin shown as percent; r-values log-transformed and scaled for visibility"
  )
