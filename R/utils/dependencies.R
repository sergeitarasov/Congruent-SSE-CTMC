# library(rphenoscate)
# library(corHMM)
# library(phytools)
# library(plyr)
# library(dplyr)
# library(diversitree)
# library(hisse)

# Install rphenoscate packages that contains functions to automatize hidden state expansions
# install.packages("devtools")
# devtools::install_github("phenoscape/rphenoscape")
# remotes::install_github("uyedaj/treeplyr")
# devtools::install_github("uyedaj/rphenoscate")
# library(rphenoscate)
# library(Rcpp)
# library(RcppArmadillo)
# library(hisse)
# library(corHMM)

# List of packages for session
.packages <- c(
  'corHMM',
  'phytools',
  'plyr',
  'dplyr',
  'hisse',
  'diversitree',
  'rphenoscape',
  'rphenoscate',
  'Rcpp',
  'partitions',
  'RcppArmadillo'
)

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)

