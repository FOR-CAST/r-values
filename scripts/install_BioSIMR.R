## NOTE: other than installing Java, this manual installation is not necessary no recommended
##       J4R and BioSimR are installed from GitHub as part of `renv` package management.

if (FALSE) {
  ## Setup (only need to run once on a new machine)
  Sys.setenv(R_LIBCURL_SSL_REVOKE_BEST_EFFORT = TRUE) #Needed for "several things that we do in the NRCan network"

  # Install j4r
  install.packages(
    "https://sourceforge.net/projects/repiceasource/files/latest/download",
    repos = NULL,
    type = "source"
  )

  # Note: On NRCan machines, must install Perforce OpenLogic OpnJDK 64-bit, as a replacement for Java
  Sys.setenv(
    PATH = paste(
      "C:/Program Files/OpenLogic/jdk-22.0.2.9-hotspot/bin",
      Sys.getenv("PATH"),
      sep = ";"
    )
  )
  options(java.home = "C:/Program Files/OpenLogic/jdk-22.0.2.9-hotspot")

  # Install BioSIM
  install.packages(
    "https://sourceforge.net/projects/biosimclient.mrnfforesttools.p/files/latest/download",
    repos = NULL,
    type = "source"
  )

  #To get past NRCan firewall:
  install.packages("remotes")
  remotes::install_github("RNCan/BioSimClient_R")
}
