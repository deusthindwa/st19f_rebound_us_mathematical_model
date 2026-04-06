# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission
# 01/09/2024

#-------------------------------------------------------------------------------

# load the require packages
if (!require(pacman)){
  install.packages("pacman")
  install.packages("RcmdrPlugin.KMggplot2")
}

pacman::p_load(
  char = c("tidyverse","remotes","dplyr","rio","tidyr","plyr","lubridate","reshape2","curl","patchwork",
           "deSolve","adaptivetau","data.table","scales","readr","MASS","rootSolve", "labelVector","PropCIs",
           "binom","coda","Rcpp","gmm","RcppArmadillo","devtools","lattice", "RColorBrewer", "lhs","png",
           "readxl","viridis","zoo","lattice","latex2exp","ape","cowplot","gridExtra","grid","ggpubr","here")
  )

