# validate and merge csv tables ---------------------------------------------------------------

csv_files <- outputPath |>
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
##
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

## cleanup any pre-existing files from previous runs
file.path(dataPath, "AB", c("csv", "gdb", "png")) |>
  purrr::walk(function(d) if (dir.exists(d)) unlink(d, recursive = TRUE))

file.path(outputPath, "AB", c("csv", "gdb", "png")) |>
  purrr::walk(function(d) if (dir.exists(d)) unlink(d, recursive = TRUE))

## read in each (set of) files, validate, and merge
log_file <- file.path(dataPath, "AB", "_mdb_data_clean_log.txt")
file.remove(log_file)

fix_coords <- function(df) {
  ## allow at most ~1 degree lon/lat error; adjust this as needed
  lat_thrsh <- 0.99
  lon_thrsh <- 0.99

  df |>
    ## longitude sometimes entered as e.g., 110 (i.e. degrees E) instead of -110 (degrees W)
    mutate(
      lon_dd = if_else(lon_dd > 0, -lon_dd, lon_dd)
    ) |>
    ## coords 0,0 were likely used as NA
    mutate(
      lat_dd = if_else(lat_dd == 0 & lon_dd == 0, NA_real_, lat_dd),
      lon_dd = if_else(lat_dd == 0 & lon_dd == 0, NA_real_, lon_dd)
    ) |>
    ## fix presumed mis-entered degrees (e.g., 50.12345 when all others in site are 54.vwxyz)
    group_by(siteID) |>
    mutate(
      lat_diff = abs(floor(lat_dd) - stat_mode(floor(lat_dd))),
      lat_outlier = lat_diff >= lat_thrsh,
      lon_diff = abs(ceiling(lon_dd) - stat_mode(ceiling(lon_dd))), ## lon is negative; use ceiling
      lon_outlier = lon_diff >= lon_thrsh
    ) |>
    group_modify(
      .f = function(.x, .y) {
        .x |>
          mutate(
            lat_dd = if_else(lat_outlier, stat_mode(floor(lat_dd)) + lat_dd %% 1, lat_dd),
            lon_dd = if_else(lon_outlier, stat_mode(ceiling(lon_dd)) - lon_dd %% 1, lon_dd)
          )
      }
    ) |>
    ungroup() |>
    ## drop accounting cols
    select(-lat_diff, -lat_outlier, -lon_diff, -lon_outlier)
}

abr <- dirname(dirname(csv_files)) |>
  unique() |>
  purrr::map(function(d) {
    cli::cli_h1(glue::glue("Processing directory {fs::path(d)}"))
    ## if 'mpb_site' table missing, use 'mpb_survey_info'
    d_site <- file.path(d, "site")
    survey_site_csv <- file.path(d_site, list.files(d_site, pattern = "_mpb_site[.]csv"))
    if (length(survey_site_csv) == 0) {
      survey_site_csv <- file.path(d_site, list.files(d_site, pattern = "_mpb_survey_info[.]csv"))
    }

    d_tree <- file.path(d, "tree")
    mpb_trees_csv <- file.path(d_tree, list.files(d_tree, pattern = "_mpb_trees[.]csv"))

    purrr::map2(survey_site_csv, mpb_trees_csv, function(fsite, ftrees) {
      ## check the csvs are correctly paired before joining
      stopifnot(identical(
        sub("(mpb_site|mpb_survey_info)", "mpb_trees", basename(fsite)),
        basename(ftrees)
      ))

      site_tbl_name <- fs::path_file(fs::path_ext_remove(fsite))
      trees_tbl_name <- fs::path_file(fs::path_ext_remove(ftrees))
      cli::cli_alert_info(glue::glue(
        "    joining tables: '{site_tbl_name}' and '{trees_tbl_name}'..."
      ))

      ## From the site files we only want:
      ## - btl_year;
      ## - siteID;
      ## - plot_lat_dd, plot_lon_dd (or plot_long_dd);
      ## - nbr_infested (or nbr_trees);
      ## - r_value;
      survey_site <- read.csv(fsite) |>
        mutate(infestation = as.character(infestation)) |>
        select(
          beetle_yr,
          siteID,
          matches("plot_(lat|lon|long)_dd"),
          matches("plot_(lat|lon|long)_dmd"),
          matches("nbr_infested|nbr_trees"),
          r_value
        ) |>
        rename(
          any_of(c(
            plot_lon_dd = "plot_long_dd", ## if plot_long_dd exists, rename to plot_lon_dd
            plot_lon_dmd = "plot_long_dmd", ## if plot_long_dmd exists, rename to plot_lon_dmd
            nbr_infested = "nbr_trees", ## if nbr_trees exists, rename to nbr_infested
            r_value_site = "r_value" ## if r_value exists, rename to r_value_site
          ))
        ) |>
        mutate(
          ## NAs generate spurious warnings
          plot_lat_dd2 = suppressWarnings(parzer::parse_lat(as.character(plot_lat_dmd))),
          plot_lon_dd2 = suppressWarnings(parzer::parse_lon(as.character(plot_lon_dmd)))
        ) |>
        mutate(plot_lat_dmd = NULL, plot_lon_dmd = NULL)

      ## bounds checking + validation:
      survey_site <- survey_site |>
        mutate(
          plot_lat_dd = if_else(plot_lat_dd == 0 & plot_lon_dd == 0, NA_real_, plot_lat_dd),
          plot_lon_dd = if_else(plot_lat_dd == 0 & plot_lon_dd == 0, NA_real_, plot_lon_dd),

          plot_lat_dd2 = if_else(plot_lat_dd2 == 0 & plot_lon_dd2 == 0, NA_real_, plot_lat_dd2),
          plot_lon_dd2 = if_else(plot_lat_dd2 == 0 & plot_lon_dd2 == 0, NA_real_, plot_lon_dd2)
        ) |>
        mutate(
          plot_lon_dd = if_else(plot_lon_dd > 0, -plot_lon_dd, plot_lon_dd),
          plot_lon_dd2 = if_else(plot_lon_dd2 > 0, -plot_lon_dd2, plot_lon_dd2)
        )

      if (grepl("mpb_survey_jun_11_2008", fsite)) {
        survey_site <- survey_site |>
          mutate(
            ## entered longitude -188; likely should be -118
            plot_lon_dd = if_else(
              siteID == 20 & floor(abs(plot_lon_dd)) == 188,
              -(118 + plot_lon_dd %% 1),
              plot_lon_dd
            ),
            plot_lon_dd2 = if_else(
              siteID == 20 & floor(abs(plot_lon_dd2)) == 188,
              -(118 + plot_lon_dd2 %% 1),
              plot_lon_dd2
            )
          )
      } else if (grepl("mpb_survey_may_01_2008", fsite)) {
        survey_site <- survey_site |>
          mutate(
            ## entered latitude 50; likely should be 51 to be in AB
            plot_lat_dd = if_else(
              siteID == 6 & round(plot_lat_dd) == 50,
              51 + plot_lat_dd %% 1,
              plot_lat_dd
            ),
            plot_lat_dd2 = if_else(
              siteID == 6 & round(plot_lat_dd2) == 50,
              51 + plot_lat_dd2 %% 1,
              plot_lat_dd2
            )
          )
      }

      survey_site <- survey_site |>
        mutate(
          ## use lat/lon_dd2 (from lat/lon_dmd) where lat/lon_dd are NA
          plot_lat_dd = if_else(
            is.na(plot_lat_dd) & !is.na(plot_lat_dd2),
            plot_lat_dd2,
            plot_lat_dd
          ),
          plot_lon_dd = if_else(
            is.na(plot_lon_dd) & !is.na(plot_lon_dd2),
            plot_lon_dd2,
            plot_lon_dd
          )
        ) |>
        mutate(plot_lat_dd2 = NULL, plot_lon_dd2 = NULL)

      ## correct for missing r-values
      if (!"r_value_site" %in% colnames(survey_site)) {
        survey_site <- mutate(survey_site, r_value_site = NA_real_)
      }

      survey_site <- survey_site |>
        mutate(
          r_value_site = if_else(r_value_site == -999, NA_real_, r_value_site)
        )

      ## beetle_yr should be 4 digits, +/- 1 year from the rest of the table
      if (length(unique(survey_site$beetle_yr)) > 1) {
        ## TODO: decide how to handle each case
        cat(
          c(
            paste(">> ", fsite),
            paste(
              "    multiple beetle years:",
              paste(unique(survey_site$beetle_yr), collapse = ", ")
            )
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
      if (is.na(btl_yr)) {
        ## extract year from the filename
        btl_yr <- stringr::str_extract(basename(fsite), "20[0-9][0-9]") |> as.integer()
      }

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
          starts_with("ns"),
          starts_with("ss"),
          matches("r_value")
        ) |>
        rename(
          any_of(c(
            lon_dd = "long_dd", ## if long_dd exists, rename to lon_dd
            r_value_tree = "r_value" ## if r_value exists, rename to r_value_tree
          ))
        ) |>
        mutate(
          ## bounds checking + validation:
          siteID = if_else(is.na(siteID), lead(siteID), as.integer(siteID)),
          tree_nbr = as.character(tree_nbr), ## character id which is sometimes a number;
          dbh = if_else(dbh > 130, 130, dbh), ## TODO: cap dbh at 130? or use idiosyncratic rules per BC email
          ht_pitch_tube = if_else(ht_pitch_tube > 25, ht_pitch_tube / 100, ht_pitch_tube) ## should be m, not cm
        )

      if (!any(c("lat_dd", "lon_dd") %in% colnames(mpb_trees))) {
        mpb_trees <- mutate(mpb_trees, lat_dd = NA_real_, lon_dd = NA_real_)
      }

      ## identify and fix erroneous tree lat/lon
      mpb_trees <- mpb_trees |> fix_coords()

      ## correct for missing r-values
      if (!"r_value_tree" %in% colnames(mpb_trees)) {
        mpb_trees <- mutate(mpb_trees, r_value_tree = NA_real_)
      }

      mpb_trees <- mpb_trees |>
        mutate(
          r_value_tree = if_else(r_value_tree == -999, NA_real_, r_value_tree)
        )

      mpb_site_trees <- left_join(mpb_trees, survey_site, by = "siteID")

      ## create maps of lat/lon plots to identify erroneous coords
      site_sf <- survey_site |>
        na.omit() |>
        sf::st_as_sf(
          coords = c("plot_lon_dd", "plot_lat_dd"),
          crs = sf::st_crs(4326)
        )

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
      p <- p +
        theme_bw(base_size = 20) +
        annotation_north_arrow(
          location = "bl",
          which_north = "true",
          pad_x = unit(0.25, "in"),
          pad_y = unit(0.25, "in"),
          style = north_arrow_fancy_orienteering
        ) +
        xlab("Longitude") +
        ylab("Latitude")

      out_path_csv <- file.path(outputPath, "AB", "csv") |> fs::dir_create()
      out_path_png <- file.path(outputPath, "AB", "png") |> fs::dir_create()

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

write.csv(
  abr,
  file = file.path(outputPath, "AB", "csv", "all_mpb_site_trees_cleaned.csv"),
  row.names = FALSE
)

## diagnostics / checking for NAs
identical(which(is.na(abr$plot_lat_dd)), which(is.na(abr$plot_lon_dd))) ## TRUE

abr_na_coords <- filter(abr, is.na(plot_lat_dd))
nrow(abr_na_coords) ## 99

all(is.na(abr_na_coords$lat_dd)) ## TRUE
all(is.na(abr_na_coords$lon_dd)) ## TRUE
