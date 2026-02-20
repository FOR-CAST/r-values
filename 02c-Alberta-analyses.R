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
      s(asin(sqrt(Q)), bs = "gp") +
      s(SSI_2023, bs = "gp") +
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

## Correct DBH values > 100 by assuming they are circumference
all_data_df_join_CMI$dbh <- ifelse(
  all_data_df_join_CMI$dbh > 100,
  all_data_df_join_CMI$dbh / (2 * pi),
  all_data_df_join_CMI$dbh
)

## below we doubled k' based on initial gam.check of model using default k-values;
## then doubled k-values for Psurv, PineVol, and SSI again.
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
    # s(PineVol, bs = "gp", k = 44) +
    s(Psurv, bs = "gp", k = 44) +
    s(SSI_2016, bs = "gp", k = 44), # +
  # s(Tmin, bs = "gp", k = 22),
  data = all_data_df_join_CMI,
  method = "REML",
  family = gaussian(link = "identity")
)
summary(gam_model.all)

gam.check(gam_model.all)

qq.gam(gam_model.all, pch = 20)

# dev.new()
plot(gam_model.all, scheme = 2, pages = 1, all.terms = TRUE)

# png(file.path(figPath, "gam_model_cleaned.png"), height = 1600, width = 1600, res = 300)
pdf(file.path(figPath, "gam_model_AB.pdf"), height = 8, width = 8)
par(cex = 1.4, cex.axis = 1.2, cex.lab = 1.4, cex.main = 1.6)
plot(gam_model.all, scheme = 2, pages = 1, all.terms = TRUE)
dev.off()

### -------------------------------------------------------------------------------------------

# ssi_crs <- get_SSI(dsn = ssi_gdb, year = 2023) |> sf::st_crs() ## EPSG:3400
ssi_crs <- sf::st_crs(3400)

## bring the geometry back into the sf
all_data_sf <- all_data_df_join_CMI |>
  dplyr::filter(!is.na(lon) & !is.na(lat) & !is.na(beetle_yr)) |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  sf::st_make_valid() |>
  sf::st_transform(ssi_crs)

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

### Plot r-values unprojected
plot_df <- all_data_df_join_CMI |>
  filter(!is.na(lon) & !is.na(lat) & !is.na(r) & !is.na(beetle_yr)) |>
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    r_log = log10(r + 1),
    beetle_yr = as.factor(beetle_yr)
  )

ab_outline_df <- ab_sf |>
  sf::st_transform(4326) |> ## match lon/lat CRS
  sf::st_coordinates() |>
  as.data.frame() |>
  rename(lon = X, lat = Y)

gg_r_by_year_unproj <- ggplot(plot_df, aes(x = lon, y = lat, fill = r_log)) +
  geom_point(size = 2, alpha = 0.8, stroke = 0.2, shape = 21, color = "black") +
  geom_path(
    data = ab_outline_df,
    aes(x = lon, y = lat),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.8
  ) +
  scale_x_continuous(
    breaks = pretty(plot_df$lon, n = 5), ## ~5 evenly spaced ticks
    labels = function(x) sprintf("%.0f", x) ## no decimals
  ) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    name = "r",
    breaks = c(0, 1, 2),
    labels = c("1", "10", "100")
  ) +
  facet_wrap(~beetle_yr, ncol = 7, nrow = 2) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude", title = "Spatial Distribution of r-values by Year")

ggsave(
  file.path(figPath, "map_r-values_by_year_unprojected.png"),
  gg_r_by_year_unproj,
  height = 9,
  width = 15
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

plot_df_summary <- yearly_summary |>
  pivot_longer(
    cols = c(mean_r_log, se_r_log, mean_Psurv, se_Psurv, mean_Tmin, se_Tmin),
    names_to = "metric",
    values_to = "value"
  )

Psurv.r.cor <- stats::cor(yearly_summary$mean_r_log, yearly_summary$mean_Psurv)

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
  geom_text(
    data = data.frame(x = 2012, y = 100, label = paste0("r² = ", round(Psurv.r.cor, 2))),
    aes(x = x, y = y, label = label),
    hjust = 1.1,
    vjust = 1.5,
    inherit.aes = FALSE,
    size = 4
  ) +
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

## Define scale and shift for log-transformed r
y_scale <- 25 ## adjust as needed for visual alignment
y_shift <- 0 ## adjust as needed for vertical offset

## Choose r values for axis tics — adjust if needed
r_breaks <- c(1, 2, 3)
r_labels <- as.character(r_breaks)
r_transformed <- log10(r_breaks) * y_scale + y_shift

gg_r_Psurv_ribbon <- ggplot(yearly_summary, aes(x = beetle_yr)) +
  ## Psurv line and ribbon
  geom_ribbon(
    aes(ymin = mean_Psurv - se_Psurv, ymax = mean_Psurv + se_Psurv, fill = "Psurv (%)"),
    alpha = 0.2,
    show.legend = FALSE
  ) +
  geom_line(aes(y = mean_Psurv, color = "Psurv (%)"), size = 1.2) +
  geom_point(
    aes(y = mean_Psurv, shape = "Psurv (%)"),
    size = 4,
    stroke = 1.2,
    color = "black",
    fill = "white"
  ) +
  geom_text(
    data = data.frame(x = 2013.5, y = 100, label = paste0("r = ", round(Psurv.r.cor, 2))),
    aes(x = x, y = y, label = label),
    hjust = 1.1,
    vjust = 1.5,
    inherit.aes = FALSE,
    size = 4
  ) +
  ## r line and ribbon (scaled)
  geom_ribbon(
    aes(
      ymin = (mean_r_log - se_r_log) * y_scale + y_shift,
      ymax = (mean_r_log + se_r_log) * y_scale + y_shift,
      fill = "r"
    ),
    alpha = 0.2,
    show.legend = FALSE
  ) +
  geom_line(aes(y = mean_r_log * y_scale + y_shift, color = "r"), size = 1.2) +
  geom_point(
    aes(y = mean_r_log * y_scale + y_shift, shape = "r"),
    size = 4,
    color = "black"
  ) +
  ## Vertical markers
  geom_vline(xintercept = c(2010, 2015), linetype = "dashed", color = "gray50") +

  ## Axes and legend
  scale_y_continuous(
    name = "Psurv (%)",
    limits = c(0, 100),
    sec.axis = sec_axis(
      transform = ~ . / 200,
      name = "r",
      breaks = c(0, 0.2998, 0.477),
      labels = c("1", "2", "3")
    )
  ) +
  geom_hline(yintercept = log10(c(1, 2, 3)) * y_scale, linetype = "dashed", color = "gray50") +

  scale_color_manual(
    values = c("Psurv (%)" = "black", "r" = "black"),
    breaks = c("Psurv (%)", "r")
  ) +
  scale_fill_manual(
    values = c("Psurv (%)" = "black", "r" = "black"),
    breaks = c("Psurv (%)", "r")
  ) +
  scale_shape_manual(
    values = c("Psurv (%)" = 22, "r" = 16),
    breaks = c("Psurv (%)", "r")
  ) +
  labs(x = "Beetle Year") +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )

ggsave(
  file.path(figPath, "mean_Psurv_r_over_time_new.png"),
  gg_r_Psurv_ribbon,
  width = 6,
  height = 4,
  dpi = 300,
  units = "in"
)

ggsave(
  file.path(figPath, "mean_Psurv_r_over_time_new.pdf"),
  gg_r_Psurv_ribbon,
  width = 6,
  height = 4
)

## Re-do the ribbon plot but without the ridiculous overlay scheme in ggplot

## Plot r on log scale (This is Figure 4 in the manuscript)
png(file.path(figPath, "r_Psurv_overtime.png"), height = 1800, width = 2400, res = 300)
par(mar = c(4, 5, 2, 6))
plot(
  yearly_summary$beetle_yr,
  10^yearly_summary$mean_r_log - 1,
  type = "b",
  pch = 16,
  col = "black",
  xlab = "Beetle Year",
  ylab = "r",
  log = "y",
  cex = 1.5
)
abline(h = 1, lty = 2, col = "red")
## Overlay Psurv on second axis
par(new = TRUE)
plot(
  yearly_summary$beetle_yr,
  yearly_summary$mean_Psurv,
  type = "b",
  pch = 22,
  bg = "white",
  cex = 1.5,
  col = "black",
  axes = FALSE,
  xlab = "",
  ylab = "",
  ylim = c(0, 100)
)
axis(4)
mtext("P (% survival)", side = 4, line = 3)
legend(
  2012,
  20,
  legend = c("r", "P"),
  pch = c(16, 22),
  pt.cex = 1.5,
  cex = 1,
  text.width = max(strwidth(c("r", "Psurv"), cex = 1))
)
text(2007, 95, paste0("r = ", round(Psurv.r.cor, 2)))
dev.off()

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

## next we make a boxplot of r-values binned by year
r.box <- ggplot(plot_df, aes(x = factor(beetle_yr), y = r_log)) +
  geom_boxplot(fill = "grey80", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.6, alpha = 0.6, color = "black") +
  geom_hline(yintercept = log10(2), color = "red", linetype = "dashed", size = 0.6) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1, 1.5, 2),
    labels = c("0", "3", "9", "31", "99"),
    name = "r"
  ) +
  labs(x = "Beetle Attack Year", title = "Distribution of r-values by Year") +
  theme_minimal()

ggsave(file.path(figPath, "boxplot_r_over_time.png"), r.box, height = 6, width = 9)

## violin plot
r.violin <- ggplot(plot_df, aes(x = factor(beetle_yr), y = r)) +
  geom_violin(fill = "grey80", color = "black") +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", size = 0.6) +
  scale_y_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100),
    labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20", "50", "100"),
    name = "r"
  ) +
  labs(x = "Beetle Attack Year", title = "Violin Plot of r-values by Year") +
  theme_minimal()

ggsave(file.path(figPath, "violinplot_r_over_time.png"), r.violin, height = 6, width = 9, dpi = 300)

## violin plot with overlaid box plot and median added
r.violin.box <- ggplot(plot_df, aes(x = factor(beetle_yr), y = r)) +
  geom_violin(fill = "grey80", color = "black") +
  stat_summary(fun.y = median, geom = "point", size = 2) +
  geom_boxplot(width = 0.1) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", size = 0.6) +
  scale_y_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100),
    labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20", "50", "100"),
    name = "r"
  ) +
  labs(x = "Beetle Attack Year", title = "Violin Plot of r-values by Year") +
  theme_minimal()

ggsave(
  file.path(figPath, "violinboxplot_r_over_time.png"),
  r.violin.box,
  height = 6,
  width = 9,
  dpi = 300
)

## proportion of zeroes (omitted from violin plot)
plot_df |>
  group_by(beetle_yr) |>
  summarise(prop_zero = mean(r == 0)) |>
  ggplot(aes(x = factor(beetle_yr), y = prop_zero)) +
  geom_col(fill = "firebrick", color = "black") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Beetle Attack Year",
    y = "Proportion of Zero r-values",
    title = "Zero Proportion by Year"
  ) +
  theme_minimal()

## Filter out zeroes for violin plot
plot_df_pos <- plot_df |> filter(r > 0)

## Compute proportion of zeroes per year
zero_prop_df <- plot_df |>
  group_by(beetle_yr) |>
  summarise(prop_zero = mean(r == 0))

## overlay proportion zeroes onto violin plot
r.violin.bar <- ggplot() +
  ## Overlay bars for proportion of zeroes
  geom_col(
    data = zero_prop_df,
    aes(x = factor(beetle_yr), y = prop_zero * max(plot_df_pos$r, na.rm = TRUE)),
    fill = "firebrick",
    alpha = 0.3,
    width = 0.6
  ) +
  ## Violin plot for r > 0
  geom_violin(
    data = plot_df_pos,
    aes(x = factor(beetle_yr), y = r),
    fill = "grey40",
    color = "black"
  ) +
  ## Boxplot overlay (same data, same x/y)
  geom_boxplot(
    data = plot_df_pos,
    aes(x = factor(beetle_yr), y = r),
    width = 0.1
  ) +

  ## Median point overlay
  stat_summary(
    data = plot_df_pos,
    aes(x = factor(beetle_yr), y = r),
    fun = median,
    geom = "point",
    size = 2
  ) +

  ## Log scale for r
  scale_y_log10(
    breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100),
    labels = c("0.05", "0.1", "0.2", "0.5", "1", "2", "5", "10", "20", "50", "100"),
    name = "r",
    sec.axis = dup_axis(name = "% zeroes")
  ) +
  labs(
    x = "Beetle Attack Year",
    title = ""
  ) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", size = 0.6) +
  geom_hline(yintercept = 0, color = "black", size = 0.8) + # solid y-axis
  geom_vline(xintercept = 0.5, color = "black", size = 0.8) + # solid x-axis at left edge
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(margin = margin(t = 10))
  )

## this is Figure 2 in the manuscript
ggsave(
  file.path(figPath, "violinbarplot_r_over_time2.png"),
  r.violin.bar,
  height = 6,
  width = 9,
  dpi = 300
)

## file RTC are red tree counts from Mike Undershultz Feb 22, 2023
rtc <- file.path(dataPath, "AB", "RedTreeCounts.txt") |>
  read.table(header = TRUE)

png(file.path(figPath, "RedGreenTrees.png"), width = 2400, height = 1800, res = 300)
par(mar = c(5, 4, 4, 6)) # bottom, left, top, right
plot(
  rtc$Year,
  log10(rtc$RedTrees),
  xlab = "Survey Year",
  yaxt = "n",
  ylab = "Thousands of Trees",
  type = "l",
  col = "red",
  ylim = c(3, 7),
  lwd = 2
)
axis(2, at = 3:7, labels = c("1", "10", "100", "1 000", "10 000"))
abline(v = c(2010, 2015, 2020), col = "gray")
points(rtc$Year, log10(rtc$RedTrees), pch = 19, col = "red", cex = 1.5)
lines(rtc$Year, log10(rtc$TreesControlled), lwd = 2, col = "darkgreen")
points(rtc$Year, log10(rtc$TreesControlled), pch = 15, col = "darkgreen", cex = 1.5)
legend(
  2013,
  7,
  legend = expression(
    "Red Trees Detected (" * X[t] * ")",
    "Green Trees Controlled",
    "Ratio of Change from Last Year (" * R[t] * ")"
  ),
  col = c("red", "darkgreen", "blue"),
  pch = c(19, 15, 17),
  lwd = 2
)

## compute and plot interannual change in red trees Rt=Xt/Xt-1, on second y-axis
rtc$Rt <- c(NA, rtc$RedTrees[-1] / rtc$RedTrees[-length(rtc$RedTrees)])

par(new = TRUE)
plot(
  rtc$Year,
  log10(rtc$Rt),
  type = "l",
  col = "blue",
  lwd = 2,
  axes = FALSE,
  xlab = "",
  ylab = "",
  ylim = c(-1, log10(50))
) ## covers ratio from ~0.1 to 50
abline(h = 0, lty = 2, col = "blue")
points(rtc$Year, log10(rtc$Rt), pch = 17, col = "blue", cex = 1.2)
axis(
  4,
  at = log10(c(0.1, 0.5, 1, 2, 10, 30, 50)),
  labels = c("0.1", "0.5", "1", "2", "10", "30", "50")
)
mtext(expression(R[t] == X[t] / X[t - 1]), side = 4, line = 3)
dev.off()

## same but in ggplot:
rtc$RedTrees_log <- log10(rtc$RedTrees)
rtc$TreesControlled_log <- log10(rtc$TreesControlled)
rtc$Rt <- c(NA, rtc$RedTrees[-1] / rtc$RedTrees[-length(rtc$RedTrees)])
rtc$Rt_log <- log10(rtc$Rt)

red.green.plot <- ggplot(rtc, aes(x = Year)) +
  ## Red Trees line and points
  geom_line(aes(y = RedTrees_log, color = "Red Trees"), size = 1.2) +
  geom_point(aes(y = RedTrees_log, color = "Red Trees", shape = "Red Trees"), size = 3) +

  ## Green Trees line and points
  geom_line(aes(y = TreesControlled_log, color = "Green Trees"), size = 1.2) +
  geom_point(aes(y = TreesControlled_log, color = "Green Trees", shape = "Green Trees"), size = 3) +

  ## Ratio line and points on secondary axis
  geom_line(aes(y = Rt_log, color = "Ratio"), size = 1.2) +
  geom_point(aes(y = Rt_log, color = "Ratio", shape = "Ratio"), size = 2.5) +

  ## Vertical reference lines
  geom_vline(xintercept = c(2010, 2015, 2020), color = "gray70", linetype = "dashed") +

  ## Horizontal reference line at ratio = 1 (log10 = 0)
  geom_hline(yintercept = 0, color = "blue", linetype = "dotted") +

  ## Primary y-axis (log10 trees)
  scale_y_continuous(
    name = "Thousands of Trees",
    breaks = 3:7,
    labels = c("1", "10", "100", "1 000", "10 000"),
    sec.axis = sec_axis(
      transform = ~., ## identity transform
      breaks = log10(c(0.1, 0.5, 1, 2, 10, 30, 50)),
      labels = c("0.1", "0.5", "1", "2", "10", "30", "50"),
      name = expression(R[t] == X[t] / X[t - 1])
    )
  ) +

  ## Color and shape mappings
  scale_color_manual(
    name = NULL,
    values = c("Red Trees" = "red", "Green Trees" = "darkgreen", "Ratio" = "blue")
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("Red Trees" = 19, "Green Trees" = 15, "Ratio" = 17)
  ) +

  labs(x = "Survey Year") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.right = element_text(color = "blue"),
    legend.position.inside = c(0.8, 0.9),
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.margin = margin(t = 10, r = 60, b = 10, l = 10) ## equivalent to par(mar)
  )

## plot Rt against rt: rtc$Rt vs plot_df$r_log -----------------------------------------------------

## Biological and Survey Timeline: Why r(t) Predicts R(t+2)
##
## This section explains the biological and observational logic behind our test:
## Does brood productivity in beetle year t (r(t)) predict outbreak growth two years later (R(t+2))?
##
## Beetle Year t:
## - Adult beetles attack trees in July of year t.
## - Eggs are laid in those trees and overwinter from fall of year t to spring of year t+1.
##
## Survey Year t+1:
## - In May–June of year t+1, ground crews use maps of red trees from the previous year t
##   to identify spots in need of ground surveys for detecting fresh green attack in need of control.
## - During these ground surveys in year t+1 we assess brood productivity (r(t)) by sampling trees attacked in year t.
## - These r-values reflect how successful the beetles were at producing viable offspring.
## - Brood maturing in year t+1 emerge as adults and disperse to new green trees, and begin a new attack cycle
## - this adult dispersal in t+1 changes the potential size of the area experiencing outbreak; it can shrink or grow
## - Later in summer (Aug–Nov of t+1), aerial surveys detect newly attacked fading trees (X(t+1)) — these are the "red trees".
## - These red trees detected in t+1 were attacked by beetles that emerged from brood sites laid in year t.
## - In BioSIM, this is the year reported for Psurv, not fall of t, but spring of t+1, the simulated year when the calculation finishes,
##   not the simulated calendar year when the calculation started
##
## Winter t+1 → t+2:
## - Control operations (cut & burn) target fading green trees with live larvae — often the same trees that will fade and turn red in year t+2.
##
## Survey Year t+2:
## - Aerial surveys detect red trees again (X(t+2)), showing the spatial expansion or contraction of the outbreak.
## - The outbreak growth rate is calculated as R(t+2) = X(t+2) / X(t+1).
##
## Causal Chain:
## - r(t) measures reproductive success of beetles in year t.
## - Those offspring emerge and attack new trees in year t+1.
## - The impact of those attacks is visible in the red tree count of year t+2.
## - Therefore, r(t) should be tested as a predictor of R(t+2).
##
## Implementation:
## - We align r(t) with R(t+2) by shifting (i.e., de-lagging) the Rt vector upward/forward by 2 years.
## - This ensures that each r(t) value is paired with the outbreak growth rate that follows two years later.
## - We do not store the shifted Rt in rtc, because rtc$Year is indexed by survey year, not beetle year.
## - Storing it there would conflate survey-year indexing with beetle-year causality, creating semantic ambiguity.
## - Instead, we place the doubly shifted Rt alongside plot_df$r_log, which is indexed by beetle year t.
## - This preserves biological alignment and causal clarity.
## - Final assignment:
##     Rt_log_delagged2 <- c(rtc$Rt_log[-(1:2)], NA, NA) has to be appended onto an annual aggregate of plot_df
##
## Note: In spreadsheet terms, this shift moves Rt "upward" to match earlier r(t).
## In R indexing, it's a forward shift — we're de-lagging Rt by 2 years to align with its presumed cause.
##
## This logic respects the biology of the beetle life cycle, the timing of surveys, and the reporting structure of BioSIM.
## It is critical to get this alignment correct to avoid false conclusions about causality.

rt.mean <- plot_df |>
  group_by(beetle_yr) |>
  summarize(r_log = mean(r_log, na.rm = TRUE)) |>
  arrange(beetle_yr)

## Shift Rt_log forward by 2 years
Rt_log_delagged2 <- c(rtc$Rt_log[-(1:2)], NA, NA)

## Match to beetle years in rt.mean
beetle_years <- as.character(rt.mean$beetle_yr) ## factor to character
Rt_log_delagged2_trimmed <- Rt_log_delagged2[rtc$Year %in% beetle_years]

## Final assignment
rt_aligned_df <- data.frame(
  beetle_yr = as.integer(beetle_years),
  r_log = rt.mean$r_log,
  Rt_log_delagged2 = Rt_log_delagged2_trimmed
)

## plot Rt (properly delagged twice) against rt
summary(lm(rt_aligned_df$Rt_log_delagged2 ~ rt_aligned_df$r_log))

png(file.path(figPath, "Ronr_byyear_Alberta.png"), width = 1800, height = 1800, res = 300)

##  Plot log-log with natural number ticks
plot(
  10^rt_aligned_df$r_log - 1,
  10^rt_aligned_df$Rt_log_delagged2,
  log = "xy", ## log scale for both axes
  xlab = expression(r[t]),
  ylab = expression(R[t + 2]),
  xlim = c(0.2, 2),
  xaxt = "n",
  yaxt = "n"
)

## Define natural number ticks
x_ticks <- c(0.2, 0.5, 1, 2)
y_ticks <- c(0.2, 0.5, 1, 2, 5, 10, 20)

## Add axes with natural number labels
axis(1, at = x_ticks, labels = x_ticks)
axis(2, at = y_ticks, labels = y_ticks)
segments(
  x0 = 0.191,
  y0 = 1,
  x1 = 1,
  y1 = 1, ## horizontal line
  col = "red",
  lwd = 2,
  lty = 2
)
segments(
  x0 = 1,
  y0 = 0.107,
  x1 = 1,
  y1 = 1, ## vertical line
  col = "red",
  lwd = 2,
  lty = 2
)

## Try excluding beetle year 2007, where the big jump on 2009 was immigration-driven

## Extract 2007 point
x_2007 <- 10^rt_aligned_df$r_log[rt_aligned_df$beetle_yr == 2007] - 1
y_2007 <- 10^rt_aligned_df$Rt_log_delagged2[rt_aligned_df$beetle_yr == 2007]

## Draw a red circle around it
symbols(
  x = x_2007,
  y = y_2007,
  circles = rep(0.03, length(x_2007)),
  inches = FALSE,
  add = TRUE,
  fg = "red",
  lwd = 2
)

## Create censored dataframe
censored_df <- rt_aligned_df[rt_aligned_df$beetle_yr != 2007, ]

## Refit model using named dataframe
Ronr.lm.censor <- lm(Rt_log_delagged2 ~ r_log, data = censored_df)
Ronr.lm.censor.sum <- summary(Ronr.lm.censor)

## Fit quadratic
Ronr_quad <- lm(Rt_log_delagged2 ~ r_log + I(r_log^2), data = censored_df)
Ronr_quad.sum <- summary(Ronr_quad)

text(0.7, 24, bquote(italic(r)^2 == .(format(Ronr_quad.sum$r.squared, digits = 3))))
text(0.7, 16, bquote(italic(p) == .(format(Ronr_quad.sum$coefficients[2, 4], digits = 3))))

## Generate x values in natural units
x_vals <- seq(min(10^censored_df$r_log - 1), max(10^censored_df$r_log - 1), length.out = 100)

## Transform to log scale for prediction
log_r_vals <- log10(x_vals + 1)

## Predict using clean newdata
y_pred <- predict(Ronr_quad, newdata = data.frame(r_log = log_r_vals))

## Back-transform y to natural units
y_vals <- 10^y_pred

## Overplot regression line
lines(x_vals, y_vals, col = "black", lwd = 2)
abline(a = 0, b = 1, col = "black", lwd = 2, lty = 2) ## add 1:1 theoretical expectation
dev.off()

## as above but with ggplot styling, similar to Mtn Parks figure

# --- 1. Fit the QUADRATIC model on log-log scale, excluding 2007 ---
mod_AB <- lm(Rt_log_delagged2 ~ r_log + I(r_log^2), data = subset(rt_aligned_df, beetle_yr != 2007))
mod_sum <- summary(mod_AB)
R2_val <- mod_sum$r.squared

## this model includes 2007
mod_AB_2007_incl <- lm(Rt_log_delagged2 ~ r_log + I(r_log^2), data = rt_aligned_df)

p_val <- mod_sum$fstatistic
p_val <- pf(p_val[1], p_val[2], p_val[3], lower.tail = FALSE)

# --- 2. Identify the 2007 immigration point ---
df_2007 <- rt_aligned_df |> filter(beetle_yr == 2007)

annot_text <- paste0(
  "R² = ",
  format(R2_val, digits = 3),
  "\n",
  "p = ",
  format(p_val, digits = 3)
)

## Center coordinates
x_center <- mean(range(rt_aligned_df$r_log))
y_center <- log10(15)

# --- 3. Build the ggplot ---
AB_Rvsr_plot <- ggplot(rt_aligned_df, aes(x = r_log, y = Rt_log_delagged2)) +

  ## (A) Solid black points
  geom_point(size = 3, color = "black") +

  ## (B) Year labels added below to ensure they overlay the shaded areas

  ## (C) Quadratic regression line with shaded SE interval
  geom_smooth(
    data = subset(rt_aligned_df, beetle_yr != 2007),
    aes(x = r_log, y = Rt_log_delagged2),
    method = "lm",
    formula = y ~ x + I(x^2),
    se = TRUE,
    color = "black",
    fill = "grey70",
    alpha = 0.4
  ) +

  ## (D) Highlight 2007 immigration point with red ring
  geom_point(
    data = df_2007,
    aes(x = r_log, y = Rt_log_delagged2),
    shape = 21,
    fill = NA,
    color = "red",
    stroke = 1.2,
    size = 5
  ) +

  ## (E) Natural-scale axis labels on log10 coordinates
  scale_x_continuous(
    name = expression(r[t]),
    breaks = log10(c(0.1, 0.2, 0.5, 1, 2) + 1),
    labels = c(0.1, 0.2, 0.5, 1, 2)
  ) +
  scale_y_continuous(
    name = expression(R[t + 2] == C[t + 2] / C[t + 1]),
    breaks = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)),
    labels = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)
  ) +
  ## Dashed regression line INCLUDING 2007
  geom_smooth(
    data = rt_aligned_df,
    aes(x = r_log, y = Rt_log_delagged2),
    method = "lm",
    formula = y ~ x + I(x^2),
    se = FALSE,
    color = "black",
    linetype = "dashed",
    linewidth = 0.8
  ) +

  ## (F) Red bounding box: r = 0.1–2, R = 0.1–2
  # annotate(
  #   "rect",
  #   xmin  = log10(0.1 + 1), xmax = log10(2 + 1),
  #   ymin  = log10(0.1),     ymax = log10(2),
  #   color = "red", fill = NA, linewidth = 1
  # ) +

  ## (G) R^2 annotation
  annotate(
    "text",
    x = x_center,
    y = y_center,
    hjust = 0.5,
    vjust = 0.5,
    label = annot_text,
    size = 5
  ) +

  ## (B) Year labels
  ggrepel::geom_text_repel(aes(label = beetle_yr), hjust = -0.5, vjust = -0.75, size = 4) +

  ## (H) Square frame and harmonized theme
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

AB_Rvsr_plot

## this is Figure 3 in the manuscript
ggsave(
  file.path(figPath, "AB_Rvsr.png"),
  AB_Rvsr_plot,
  height = 6,
  width = 6
)

## Replot but with Mtn Parks contrast as overlay
## This is Figure 6 in the manuscript

## Rebuild the MtnParks model on log R
## You must have already run the Jasper script to get this data structure
MtnParks_plot_df <- read.csv(file.path(outputPath, "MtnParks_plot_df.csv"))
MtnParks_mod <- lm(Rt_log ~ r_log, data = MtnParks_plot_df)
MtnParks_mod_sum <- summary(MtnParks_mod)

mp.r2 <- format(round(MtnParks_mod_sum$r.squared, 2), nsmall = 2)
mp.pval <- formatC(MtnParks_mod_sum$coefficients[2, "Pr(>|t|)"], format = "f", digits = 5)
annot_text_Mtn <- paste0("R² = ", mp.r2, "\np = ", mp.pval)

## sanity test plot
sanity.test.plot <- ggplot(MtnParks_plot_df, aes(x = r_log, y = Rt_log)) +
  geom_point(size = 3, color = "#6A1B1A") +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "#6A1B1A",
    fill = "#6A1B1A",
    alpha = 0.2
  ) +

  ## --- X AXIS: log scale, natural labels (IDENTICAL TO ALBERTA) ---
  scale_x_continuous(
    breaks = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)),
    labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20"),
    name = "rₜ (natural scale)"
  ) +

  ## --- Y AXIS: log scale, natural labels (IDENTICAL TO ALBERTA) ---
  scale_y_continuous(
    breaks = log10(c(0.1, 0.2, 0.5, 1, 2, 5)),
    labels = c("0.1", "0.2", "0.5", "1", "2", "5"),
    name = "Rₜ₊₂ (natural scale)"
  ) +

  ggrepel::geom_text_repel(aes(label = beetle_yr), vjust = -1, size = 4) +

  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

ABvsMtnParks_Rvsr_plot <- ggplot(rt_aligned_df, aes(x = r_log, y = Rt_log_delagged2)) +
  ## (A) Solid black points
  geom_point(size = 3, color = "black") +

  ## (B) Year labels added below to ensure they overlay the shaded areas

  ## (C) Quadratic regression line with shaded SE interval
  geom_smooth(
    data = subset(rt_aligned_df, beetle_yr != 2007),
    aes(x = r_log, y = Rt_log_delagged2),
    method = "lm",
    formula = y ~ x + I(x^2),
    se = TRUE,
    color = "black",
    fill = "grey70",
    alpha = 0.4
  ) +

  ## (D) Highlight 2007 immigration point with red ring
  geom_point(
    data = df_2007,
    aes(x = r_log, y = Rt_log_delagged2),
    shape = 21,
    fill = NA,
    color = "red",
    stroke = 1.2,
    size = 5
  ) +

  ## (E) Natural-scale axis labels on log10 coordinates
  scale_x_continuous(
    name = expression(r[t]),
    breaks = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20) + 1),
    labels = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20),
    limits = c(log10(0.2 + 1), log10(16 + 1))
  ) +
  scale_y_continuous(
    name = expression(R[t + 2] == C[t + 2] / C[t + 1]),
    breaks = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)),
    labels = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)
  ) +

  ## (G) R^2 annotation
  annotate(
    "text",
    x = x_center,
    y = y_center,
    hjust = 0.5,
    vjust = 0.5,
    label = annot_text,
    size = 5
  ) +

  ## (B) Year labels
  ggrepel::geom_text_repel(aes(label = beetle_yr), vjust = -1, size = 4) +

  ## (H) Square frame and harmonized theme
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

## combine the MtnParks and AB datasets to simplify plotting with ggplot2
ABvsMtnParks_Rvsr_df <- bind_rows(
  MtnParks_plot_df |> select(beetle_yr, r_log, Rt_log) |> mutate(source = "MtnParks"),
  rt_aligned_df |> rename(Rt_log = Rt_log_delagged2) |> mutate(source = "AB")
) |>
  group_by(source)

ABvsMtnParks_Rvsr_plot.orig

ABvsMtnParks_Rvsr_plot <- ggplot(ABvsMtnParks_Rvsr_df, aes(x = r_log, y = Rt_log)) +
  ## plot all points; Mountain parks in maroon, AB in black/grey
  geom_point(aes(col = source), size = 3) +
  scale_color_manual(values = c(MtnParks = "#6A1B1A", AB = "black")) +

  ## add circle around the 2007 data point in the top left corner
  geom_point(
    data = subset(ABvsMtnParks_Rvsr_df, source == "AB" & beetle_yr == 2007),
    aes(x = r_log, y = Rt_log),
    shape = 21,
    fill = NA,
    color = "red",
    stroke = 1.2,
    size = 5
  ) +

  ## add shaded regions and trendlines
  geom_smooth(
    data = subset(ABvsMtnParks_Rvsr_df, source == "MtnParks"),
    aes(col = source, fill = source),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    alpha = 0.2
  ) +
  geom_smooth(
    data = subset(ABvsMtnParks_Rvsr_df, source == "AB" & beetle_yr != 2007),
    aes(col = source, fill = source),
    method = "lm",
    formula = y ~ x + I(x^2),
    se = TRUE,
    alpha = 0.4
  ) +
  scale_fill_manual(values = c(MtnParks = "#6A1B1A", AB = "grey20")) +

  ## set axes EXACTLY as in your Alberta block
  scale_x_continuous(
    name = expression(r[t]),
    breaks = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20) + 1),
    labels = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20),
    limits = c(log10(0.2 + 1), log10(20 + 1))
  ) +
  scale_y_continuous(
    name = expression(R[t + 2] == C[t + 2] / C[t + 1]),
    breaks = log10(c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)),
    labels = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20)
  ) +

  ## R² annotations
  annotate(
    "text",
    x = x_center,
    y = y_center,
    hjust = 0.5,
    vjust = 0.5,
    label = annot_text, ## Alberta
    size = 5
  ) +
  annotate(
    "text",
    x = log10(5),
    y = log10(0.08),
    hjust = 0.5,
    vjust = 0.5,
    label = annot_text_Mtn, ## Mtn Parks
    size = 5,
    col = "#6A1B1A"
  ) +

  ## add years as labels for points
  ggrepel::geom_text_repel(
    data = ABvsMtnParks_Rvsr_df,
    aes(x = r_log, y = Rt_log, label = beetle_yr, color = source),
    min.segment.length = 0,
    hjust = -0.75,
    vjust = -1,
    size = 4
  ) +

  ## theme (unchanged)
  theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

## This is Figure 6 in the manuscript
ggsave(
  file.path(figPath, "MtnParksvsAB_Rvsr.png"),
  ABvsMtnParks_Rvsr_plot,
  height = 6,
  width = 6,
  dpi = 300
)
