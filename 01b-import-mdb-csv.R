# packages ------------------------------------------------------------------------------------

# library(archive)
library(dplyr)
# library(fs)
# library(geodata)
library(ggplot2)
library(ggspatial)
# library(googledrive)
# library(purrr)
# library(sf)

# setup ---------------------------------------------------------------------------------------

## paths
dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
figPath <- "figures" |> fs::dir_create()
outputPath <- "outputs" |> fs::dir_create()

# geospatial objects for plotting -------------------------------------------------------------

ab_sf <- geodata::gadm("CAN", level = 1, path = dataPath) |>
  sf::st_as_sf() |>
  filter(NAME_1 == "Alberta") |>
  sf::st_geometry()

# validate and merge csv tables ---------------------------------------------------------------

csv_zip <- file.path(dataPath, "extracted_mdb_tables.zip")

if (!file.exists(csv_zip)) {
  googledrive::as_id("16OK48u6g5-JdWvEFhIJ0vMxvcK-k6z5l") |>
    googledrive::drive_download(path = csv_zip)
  archive::archive_extract(csv_zip, dataPath)
}

csv_files <- dataPath |>
  list.files(pattern = "_mpb_(site|trees)[.]csv", full.names = TRUE, recursive = TRUE) |>
  fs::path_rel()

#' Calculate the statistical mode
#'
#' Based on <https://stackoverflow.com/a/8189441/1380598>
#'
stat_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

## NOTES:
##
## 1. The names of the geospatials differ between site & tree files (site files pre-pend _dd with "plot").
##    When we join the tables we want to keep both sets of geospatials, site and tree.
##    We are keeping both sets of geospatials is in case siteIDs are erroneous, missing or don't match.
##    We can always associate site with tree by geospatial proximity if siteID ever fails.
##    Also, if lat/lon in one file are erroneous (which they are in rare cases), they can be guessed at from the other.
##    We do not want to lose a single r_value because of simple metadata snafu.

## 2. Both files contain 'r_value' under the same name. We will not analyze either,
##    but we want to examine the consequences of computing r_value different ways, including their
##    two ways, using whatever method they used. We what to know if method matters.
##    We also want to know if scale matters. So we will want to name these pre-computed variables
##    "r_value_site" and "r_value_tree", as distinguished from the "r_value" that we
##    will compute from the "ns1" etc. data.
##    In our computation, we will test what happens when we alter the denominator policy.
##    Dividing by zero holes gives an Inf, which can't be used in a model.
##    But zero holes happens more often than you'd like, because the disks are only 4" circles, not 12" squares.
##    Adding 1 to every denominator biases the r-value low, however,
##    (a) the bias is small (and measurable); and (b) we don't care about absolute values.
##    We care about spatial variation in pattern. So adding 1 to all denominators avoids having to
##    chuck hard-earned count data because of a spurious Inf.
##    We do not want to lose a single r_value because of simple sampling snafu.

## read in each (set of) files, validate, and merge
log_file <- file.path(dataPath, "AB", "_mdb_data_clean_log.txt")
file.remove(log_file)

all_data <- dirname(dirname(csv_files)) |>
  unique() |>
  purrr::map(function(d) {
    ## if mpb_site table missing, use mpb_survey_info
    d_site <- file.path(d, "site")
    survey_site_csv <- file.path(d_site, list.files(d_site, pattern = "_mpb_site[.]csv"))
    if (length(survey_site_csv) == 0) {
      survey_site_csv <- file.path(d_site, list.files(d_site, pattern = "_mpb_survey_info[.]csv"))
    }

    d_tree <- file.path(d, "tree")
    mpb_trees_csv <- file.path(d_tree, list.files(d_tree, pattern = "_mpb_trees[.]csv"))

    purrr::map2(survey_site_csv, mpb_trees_csv, function(fsite, ftrees) {
      ## check the csvs are correctly paired before joining
      stopifnot(identical(sub("(mpb_site|mpb_survey_info)", "mpb_trees", basename(fsite)), basename(ftrees)))

      ## From the site files we only want:
      ## - btl_year;
      ## - siteID;
      ## - plot_lat_dd, plot_lon_dd (or plot_long_dd);
      ## - nbr_infest (or nbr_trees);
      ## - r_value;
      survey_site <- read.csv(fsite) |>
        mutate(infestation = as.character(infestation)) |>
        select(
          beetle_yr, siteID,
          matches("plot_(lat|lon|long)_dd"),
          matches("nbr_infest|nbr_trees"),
          r_value
        ) |>
        rename(
          any_of(c(
            plot_lon_dd = "plot_long_dd", ## if plot_long_dd exists, rename to plot_lon_dd
            nbr_trees = "nbr_infest",     ## if nbr_infest exists, rename to nbr_trees
            r_value_site = "r_value"      ## if r_value exists, rename to r_value_site
          ))
        )

      ## bounds checking + validation:
      if (!"r_value_site" %in% colnames(survey_site)) {
        survey_site <- mutate(survey_site, r_value_site = NA_real_)
      }

      ## beetle_yr should be 4 digits, +/- 1 year from the rest of the table
      if (length(unique(survey_site$beetle_yr)) > 1) {
        ## TODO: decide how to handle each case
        cat(
          c(
            paste(">> ", fsite),
            paste("    multiple beetle years:", paste(unique(survey_site$beetle_yr), collapse = ", "))
          ),
          file = log_file,
          append = TRUE,
          sep = "\n"
        )
        if (any(abs(diff(survey_site$beetle_yr)) > 1)) {
          cat(
            c(
              paste("    setting all beetle years to:", stat_mode(survey_site$beetle_yr))
            ),
            file = log_file,
            append = TRUE,
            sep = "\n"
          )
          survey_site$beetle_yr <- stat_mode(survey_site$beetle_yr)
        }
      }

      btl_yr <- stat_mode(survey_site$beetle_yr)

      ## From the trees files we only want:
      ## - siteID;
      ## - tree_nbr;
      ## - lat_dd, lon_dd (or long_dd);
      ## - dbh;
      ## - ht_pitch_tube;
      ## - all the survival numbers (ns1, ss1, etc.);
      ## - r_value;
      mpb_trees <- read.csv(ftrees) |>
        select(
          siteID,
          tree_nbr,
          matches("(lat|lon|long)_dd"),
          dbh,
          ht_pitch_tube,
          starts_with("ns"), starts_with("ss"),
          matches("r_value")
        ) |>
        rename(
          any_of(c(
            lon_dd = "long_dd",      ## if long_dd exists, rename to lon_dd
            r_value_tree = "r_value" ## if r_value exists, rename to r_value_tree
          ))
        ) |>
        mutate(
          ## bounds checking + validation:
          siteID = as.integer(siteID),
          tree_nbr = as.character(tree_nbr), ## character id which is sometimes a number;
          ht_pitch_tube = case_when(ht_pitch_tube > 25 ~ ht_pitch_tube / 100) ## should be m, not cm
        )

      if (!any(c("lat_dd", "lon_dd") %in% colnames(mpb_trees))) {
        mpb_trees <- mutate(mpb_trees, lat_dd = NA_real_, lon_dd = NA_real_)
      }

      if (!"r_value_tree" %in% colnames(mpb_trees)) {
        mpb_trees <- mutate(mpb_trees, r_value_tree = NA_real_)
      }

      mpb_site_trees <- left_join(mpb_trees, survey_site, by = "siteID")

      ## create maps of lat/lon plots to identify erroneous coords
      site_sf <- sf::st_as_sf(survey_site, coords = c("plot_lon_dd", "plot_lat_dd"), crs = sf::st_crs(4326))

      if (all(is.na(mpb_trees$lat_dd)) || all(is.na(mpb_trees$lon_dd))) {
        p <- ggplot() +
          geom_sf(data = ab_sf) +
          geom_sf(data = site_sf, aes(col = as.factor(siteID)))
      } else {
        trees_sf <- sf::st_as_sf(mpb_trees, coords = c("lon_dd", "lat_dd"), crs = sf::st_crs(4326))

        p <- ggplot() +
          geom_sf(data = ab_sf) +
          geom_sf(data = site_sf, aes(col = as.factor(siteID))) +
          geom_sf(data = trees_sf, aes(col = as.factor(siteID)))
      }
      p <- p  +
        theme_bw(base_size = 20) +
        annotation_north_arrow(
          location = "bl", which_north = "true",
          pad_x = unit(0.25, "in"), pad_y = unit(0.25, "in"),
          style = north_arrow_fancy_orienteering
        ) +
        xlab("Longitude") +
        ylab("Latitude")

      out_path_csv <- file.path(dataPath, "AB", "csv") |> fs::dir_create()
      out_path_png <- file.path(dataPath, "AB", "png") |> fs::dir_create()

      fname_csv <- basename(fsite) |>
        sub("(mpb_site|mpb_survey_info)", "mpb_site_trees_cleaned", x = _)
      fname_csv <- paste0(btl_yr, "_", fname_csv) |> gsub(" ", "_", x = _)

      write.csv(mpb_site_trees, file = file.path(out_path_csv, fname_csv), row.names = FALSE)

      fname_png <- paste0(tools::file_path_sans_ext(fname_csv), ".png")
      ggsave(plot = p, filename = file.path(out_path_png, fname_png))

      ## capture and log ggplot2 warnings
      cat(paste(">> plotting coords for", fsite), file = log_file, append = TRUE, sep = "\n")
      warnings(file = log_file, sep = "\n", append = TRUE)

      mpb_site_trees
    }) |>
      purrr::list_rbind()
}) |>
  purrr::list_rbind()

write.csv(all_data, file = file.path(dataPath, "AB", "csv", "all_mpb_site_trees_cleaned.csv"), row.names = FALSE)
