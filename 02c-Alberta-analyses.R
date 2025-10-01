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
  dplyr::filter(!is.na(lon) & !is.na(lat) & !is.na(beetle_yr)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_make_valid() |>
  st_transform(ssi_crs)

gg_r_by_year <- ggplot(all_data_sf) +
  geom_sf(aes(color = log10(r + 1)), linewidth = 1.1, alpha = 0.7) +
  geom_sf(data = ab_sf, fill = NA) +
  facet_wrap(~beetle_yr, ncol = 7, nrow = 2) +
  scale_color_viridis_c(option = "plasma", name = "log₁₀(r + 1)") +
  theme_minimal() +
  labs(
    title = "Spatial Distribution of r-values by Year",
    subtitle = "Log-scaled to handle extreme values",
    caption = "Each point represents a tree-level estimate"
  )

ggsave(file.path(figPath, "map_r-values_by_year.png"), gg_r_by_year, height = 16, width = 16)

### Plot r-values in base R, unprojected
plot_df <- all_data_df_join_CMI |>
  filter(!is.na(lon) & !is.na(lat) & !is.na(r) & !is.na(beetle_yr)) |>
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    r_log = log10(r + 1),
    beetle_yr = as.factor(beetle_yr)
  )

ab_outline_df <- ab_sf |>
  st_transform(4326) |> ## match lon/lat CRS
  st_coordinates() |>
  as.data.frame() |>
  rename(lon = X, lat = Y)

gg_r_by_year_unproj <- ggplot(plot_df, aes(x = lon, y = lat, fill = r_log)) +
  geom_point(size = 2, alpha = 0.8, stroke = 0.2, shape = 21, color = "black") +
<<<<<<< Updated upstream
  geom_path(
    data = ab_outline_df,
    aes(x = lon, y = lat),
    inherit.aes = FALSE,
    color = "black",
    size = 0.8
  ) +
  scale_fill_gradient(low = "white", high = "black", name = "log₁₀(r + 1)") +
  facet_wrap(~beetle_yr, ncol = 7, nrow = 2) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude", title = "Spatial Distribution of r-values by Year")

ggsave(
  file.path(figPath, "map_r-values_by_year_unprojected.png"),
  gg_r_by_year_unproj,
  height = 8,
  width = 16
)
=======
  geom_path(data = ab_outline_df, aes(x = lon, y = lat), inherit.aes = FALSE, color = "black", size = 0.8) +
  scale_x_continuous(
    breaks = pretty(plot_df$lon, n = 5),  # ~5 evenly spaced ticks
    labels = function(x) sprintf("%.0f", x)  # no decimals
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    name = "r",
    breaks = c(0, 1, 2),
    labels = c("1", "10", "100")
  ) +
  facet_wrap(~ beetle_yr, ncol = 7, nrow = 2) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude", title = "Spatial Distribution of r-values by Year")

ggsave(file.path(figPath, "map_r-values_by_year_unprojected.png"), gg_r_by_year_unproj, height = 9, width = 15)
>>>>>>> Stashed changes

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

y_scale <- 200
y_shift <- 0
gg_r_Psurv_ribbon <- ggplot(yearly_summary, aes(x = beetle_yr)) +
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
      ymin = (mean_r_log - se_r_log) * y_scale + y_shift,
      ymax = (mean_r_log + se_r_log) * y_scale + y_shift,
      fill = "log₁₀(r + 1)"
    ),
    alpha = 0.2
  ) +
  geom_line(aes(y = mean_r_log * y_scale, color = "log₁₀(r + 1)"), size = 1.2) +
  geom_point(aes(y = mean_r_log * y_scale, color = "log₁₀(r + 1)"), size = 2) +

  ## axes and legend
  scale_y_continuous(
    name = "Psurv (%)",
    limits = c(0, 100),
    sec.axis = sec_axis(~ (. - y_shift) / y_scale, name = "log₁₀(r + 1)")
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

ggsave(file.path(figPath, "mean_Psurv_r_over_time.png"), gg_r_Psurv_ribbon)

y_scale <- 75
y_shift <- -50
gg_r_Tmin_ribbon <- ggplot(yearly_summary, aes(x = beetle_yr)) +
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
    limits = c(-45, -15),
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

ggsave(file.path(figPath, "mean_Tmin_r_over_time.png"), gg_r_Tmin_ribbon)

#Next we make a boxplot of r-values binned by year
r.box<-ggplot(plot_df, aes(x = factor(beetle_yr), y = r_log)) +
  geom_boxplot(fill = "grey80", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.6, alpha = 0.6, color = "black") +
  geom_hline(yintercept = log10(2), color = "red", linetype = "dashed", size = 0.6) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "2", "9", "31", "99"),
    name = "r"
  ) +
  labs(x = "Year", title = "Distribution of r-values by Year") +
  theme_minimal()

ggsave(file.path(figPath, "boxplot_r_over_time.png"), r.box, height = 6, width = 9)

# File RTC are red tree counts from Mike Undershultz Feb 22, 2023
rtc <- read.table("data/ab/RedTreeCounts.txt", header = T)

par(mar = c(5, 4, 4, 6))  # bottom, left, top, right
plot(rtc$Year, log10(rtc$RedTrees), xlab = "Survey Year",yaxt = "n",
     ylab = "Thousands of Trees", type = "l", col = "red", ylim = c(3, 7),lwd=2)
axis(2, at = 3:7, labels = c("1", "10", "100", "1 000", "10 000"))
abline(v=c(2010,2015,2020),col="gray")
points(rtc$Year, log10(rtc$RedTrees), pch = 19, col = "red", cex = 1.5)
lines(rtc$Year, log10(rtc$TreesControlled), lwd=2, col = "darkgreen")
points(rtc$Year, log10(rtc$TreesControlled), pch = 15, col = "darkgreen", cex = 1.5)
legend(2014, 7,
       legend = expression("Red Trees Detected ("*X[t]*")",
                           "Green Trees Controlled",
                           "Ratio of Change from Last Year ("*R[t]*")"),
       col = c("red", "darkgreen", "blue"),
       pch = c(19, 15, 17),
       lwd = 2)

#compute and plot interannual change in red trees Rt=Xt/Xt-1, on second y-axis
rtc$Rt <- c(NA, rtc$RedTrees[-1] / rtc$RedTrees[-length(rtc$RedTrees)])

par(new = TRUE)
plot(rtc$Year, log10(rtc$Rt),
     type = "l", col = "blue", lwd = 2,
     axes = FALSE, xlab = "", ylab = "",
     ylim = c(-1, log10(50)))  # covers ratio from ~0.1 to 50
abline(h=0,lty=2,col="blue")
points(rtc$Year, log10(rtc$Rt), pch = 17, col = "blue", cex = 1.2)
axis(4, at = log10(c(0.1, 0.5, 1, 2, 10, 30, 50)),
     labels = c("0.1", "0.5", "1", "2", "10", "30", "50"))
mtext(expression(R[t] == X[t] / X[t-1]), side = 4, line = 3)
