# packages ------------------------------------------------------------------------------------

# library(archive)
library(dplyr)
# library(fs)
# library(googledrive)
# library(purrr)

# setup ---------------------------------------------------------------------------------------

## paths
dataPath <- normalizePath("./data", mustWork = FALSE) |> fs::dir_create()
figPath <- "figures" |> fs::dir_create()
outputPath <- "outputs" |> fs::dir_create()

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

#Notes

#From the trees files we only want:
#-siteID
#-lat_dd, lon_dd
#-dbh
#-ht_pitch_tube
#-all the survival numbers (ns1 etc)
#-r_value
#We don't care about crew, access, or any of the district information.

#From the site files we only want:
#-siteID
#-plot_lat_dd
#-plot_lon_dd
#-nbr_infest (or nbr_trees)
#-r_value
#We don't care about crew, access, or any of the district information.

#1. Note that the names of the geospatials differ between site & tree files. (site files pre-pend _dd with "plot".)
#   Good. When we join the tables we want to keep both sets of geospatials, site and tree.
#   So we can use the names as they exist. The reason we are
#   keeping both sets of spatials is in case siteIDs are erroneous, missing or don't match.
#   We can always associate site with tree by geospatial proximity if siteID ever fails.
#   Also, if lat/lons in one file are erroneous, they can be guessed at from the other.
#   Why are we doing this? Because we do not want to lose a single r_value because of simple metadata snafu.

#2. Note that r_value exists in both files under the same name. We will not analyze either,
#   but we want to examine the consequences of computing r_value different ways, including their
#   two ways, using whatever method they used. We what to know if method matters.
#   We also want to know if scale matters. So we will want to name these pre-computed variables
#   "r_value_site" and "r_value_tree", as distinguished from the "r_value" that we
#   will compute from the "ns1" etc. data. In our computation,
#   we will test what happens when we alter the denominator policy. Dividing by zero holes gives an Inf, which can't be used in a model.
#   But zero holes happens more often than you'd like, because the disks are only 4" circles, not 12" squares. Adding 1 to every
#   every denominator biases the r-value low ... but (a) the bias is small (and measurable), and (b) we don't care about absolute values
#   in our model. We care about spatial variation in pattern.


## read in each (set of) files, validate, and merge
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
browser()
      survey_site <- read.csv(fsite) |>
        mutate(infestation = as.character(infestation))

      btl_yr <- unique(survey_site$beetle_yr)

      mpb_trees <- read.csv(ftrees) |>
        mutate(infestation = as.character(infestation))

      mpb_site_trees <- left_join(mpb_trees, survey_site, by = "siteID") |>
        rename(infestation = infestation.y, site_nbr = site_nbr.y) |>
        mutate(infestation.x = NULL, site_nbr.x = NULL) |>
        select(
          beetle_yr
          ## TODO: only take the columns we care about
        ) |>
        mutate(
          ## type checking
          access_lat_dmd = as.character(access_lat_dmd), ## TODO: verify format - why *lat_dmd have 2 values? REPLY: These are degrees minutes decimal (dmd). We want decimal degrees (dd).
          access_long_dmd = as.character(access_long_dmd), ## TODO: verify format - why *lond_dmd have 2 values?
          corp_area = as.character(corp_area),
          project = as.character(project),
          project_mgr = as.character(project_mgr),
          sample_tree_nbr = as.integer(sample_tree_nbr),
          site_nbr = as.integer(site_nbr),
          surv_date = as.character(surv_date), ## TODO: use date format?
          tree_nbr = as.integer(tree_nbr)
          ## TODO: perform additional validation
        )

      fname_new <- basename(fsite) |>
        sub("(mpb_site|mpb_survey_info)", "mpb_site_trees_cleaned", x = _)
      fname_new <- paste0(btl_yr, "_", fname_new)

      file.path(dataPath, "AB", "csv") |>
        fs::dir_create() |>
        file.path(fname_new) |>
        write.csv(mpb_site_trees, file = _, row.names = FALSE)

      mpb_site_trees
    }) |>
      purrr::list_rbind()
}) |>
  purrr::list_rbind()

