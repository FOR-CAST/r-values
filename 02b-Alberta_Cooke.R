# setup ---------------------------------------------------------------------------------------

## paths
dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
figPath <- "figures" |> fs::dir_create()
outputPath <- "outputs" |> fs::dir_create()

abr<-read.csv(file.path(outputPath, "AB", "csv", "all_mpb_site_trees_cleaned.csv"), header = T)

#data check
win.graph()
par(mfrow=c(2,4))
plot(abr$lon_dd,abr$lat_dd,col="red")
plot(abr$plot_lon_dd,abr$plot_lat_dd,col="blue")
hist(abr$beetle_yr,breaks=c(2006:2019))
hist(log10(abr$nbr_infested+1))
hist(abr$dbh[abr$dbh<100])
hist(abr$ht_pitch_tube)
hist(log10(abr$r_value_site)+1)
hist(log10(abr$r_value_tree)+1)

h <- hist(abr$nbr_infested, plot = FALSE)
plot(h, yaxt = "n", main = "Histogram with Log-Scaled Y-Axis", xlab = "nbr_infested")
axis(2, at = pretty(log10(h$counts)), labels = 10^pretty(log10(h$counts)))

table(abr$beetle_yr)

his(all_abrdata$dbh[!is.na(abr$dbh)])

sum(abr$r_value_site[!is.na(abr$r_value_site)]==-999)
sum(abr$r_value_tree[!is.na(abr$r_value_tree)]==-999)

abr$r_value_site[all_data$r_value_site[!is.na(all_data$r_value_site)]==-999]

abr$r_value_site[all_data$r_value_site[!is.na(all_data$r_value_site)]==-999]<-NA

abr$r_value_site==-999

str(abr)

table (abr$nbr_infested)

library(dplyr)

abr$live <- rowSums(abr[, c(
  "ns1_eggs_live",
  "ns1_larvae_live",
  "ns1_pupae_live",
  "ns1_teneral_adults_live",
  "ns2_eggs_live",
  "ns2_larvae_live",
  "ns2_pupae_live",
  "ns2_teneral_adults_live",
  "ss1_eggs_live",
  "ss1_larvae_live",
  "ss1_pupae_live",
  "ss1_teneral_adults_live",
  "ss2_eggs_live",
  "ss2_larvae_live",
  "ss2_pupae_live",
  "ss2_teneral_adults_live"
  )], na.rm = TRUE)

win.graph()
hist(log10(abr$live+1))

abr$holes <- rowSums(abr[, c(
  "ns1_holes",
  "ns2_holes",
  "ss1_holes",
  "ss2_holes"
)], na.rm = TRUE)
win.graph()
hist(log10(abr$holes+1))

abr$r<-abr$live/(abr$holes+1) #adding 1 in denominator avoids a meaningless division by zero "error", with a relatively small cost in basing the r-value low
win.graph()
hist(log10(abr$r+1)) #lots of zero r-values: avoid log of zeroes

# Create histogram object without plotting
hist_obj <- hist(log10(abr$r+1), plot = FALSE)

mids <- hist_obj$mids

# Assign colors to bins based on average raw r-values in each bin (using BC FIDS heuristic)
bin_colors <- sapply(1:length(hist_obj$breaks[-1]), function(i) {
  # Get indices of values in this bin
  in_bin <- log10(abr$r+1) >= hist_obj$breaks[i] & log10(abr$r+1) < hist_obj$breaks[i + 1]
  avg_r <- mean(abr$r[in_bin], na.rm = TRUE)

  #BC MOF FIDS heuristic
  if (is.na(avg_r)) {
    return("white")  # default if no data in bin
  } else if (avg_r < 1) {
    return("white")
  } else if (avg_r < 2) {
    return("grey")
  } else {
    return("black")
  }
})

# Plot manually using rect()
plot(hist_obj, col = NA, border = NA, main = "Histogram of log10(r + 1)", xlab = "log10(r + 1)")
for (i in 1:length(hist_obj$counts)) {
  rect(xleft = hist_obj$breaks[i],
       xright = hist_obj$breaks[i + 1],
       ybottom = 0,
       ytop = hist_obj$counts[i],
       col = bin_colors[i],
       border = "black")
}

legend("topright",
       legend = c("r < 1: Declining", "1 ≤ r < 2: Stable", "r ≥ 2: Growing"),
       fill = c("white", "grey", "black"),
       border = "black",
       title = "BC FIDS Population Phase")

abr.early<-abr[(abr$beetle_yr<2016),]
abr.late<-abr[(abr$beetle_yr>2015),]

sapply(abr.early[, c("r", "beetle_yr", "plot_lat_dd", "plot_lon_dd", "nbr_infested", "dbh", "ht_pitch_tube")], function(x) sum(is.na(x)))
sapply(abr.late[, c("r", "beetle_yr", "plot_lat_dd", "plot_lon_dd", "nbr_infested", "dbh", "ht_pitch_tube")], function(x) sum(is.na(x)))

abr.lm.early <- lm(log10(r + 1) ~ beetle_yr + plot_lat_dd + plot_lon_dd +
                     log10(nbr_infested + 1) * dbh * ht_pitch_tube,
                   data = abr.early,
                   na.action = na.omit)
summary(abr.lm.early)
abr.lm.late <- lm(log10(r + 1) ~ beetle_yr + plot_lat_dd + plot_lon_dd +
                     log10(nbr_infested + 1) * dbh * ht_pitch_tube,
                   data = abr.late,
                   na.action = na.omit)
summary(abr.lm.late)

