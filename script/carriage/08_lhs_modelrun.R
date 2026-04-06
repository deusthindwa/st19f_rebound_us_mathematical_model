# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission dynamics and impact of vaccination policies
# 01/09/2024

#====================================================================

#load R packages
source(here::here("script", "00_main.R"))

#load data to fit model to
source(here::here("script", "01_preprocessData.R"))

#set up model including births, contacts, population
source(here::here("script", "02_setUp.R"))

#set up fixed parameter values
source(here::here("script", "03_modelParKnown.R"))

#run model function
source(here::here("script", "04_transModel.R"))

# model fit using MAP (disabled)
source(here::here("script", "05_modelfit_MAP.R"))

#call latin hypercube samples
source(here::here("script", "06_lhs_paraspace.R"))

#create a list to hold parameter estimates
saved_carr_dynamics <- list()

#number of model runs
N_run <- 1

for( nn in 1:len_each_run){
  
  #parameter sets to run the dynamic model
  parmAdded <- list(baseV  = c(lhs.sample.mt[nn,1], lhs.sample.mt[nn,2], lhs.sample.mt[nn,3], lhs.sample.mt[nn,4]),
                    baseF  = c(lhs.sample.mt[nn,5], lhs.sample.mt[nn,6], lhs.sample.mt[nn,7], lhs.sample.mt[nn,8]),
                    baseN  = c(lhs.sample.mt[nn,9], lhs.sample.mt[nn,10], lhs.sample.mt[nn,11], lhs.sample.mt[nn,12])
                    # compV  = c(lhs.sample.mt[nn,5], lhs.sample.mt[nn,6], lhs.sample.mt[nn,7], lhs.sample.mt[nn,8]),
                    # compF  = c(lhs.sample.mt[nn,9], lhs.sample.mt[nn,10], lhs.sample.mt[nn,11], lhs.sample.mt[nn,12]),
                    # compN  = c(lhs.sample.mt[nn,13], lhs.sample.mt[nn,14], lhs.sample.mt[nn,15], lhs.sample.mt[nn,16])
                    )
  
  #add 'parms' without parameters to estimate and add the parm_added with parameters to estimate
  #my_parmset has both fixed and estimated' parameters
  my_parmset <- c(parms, parmAdded)
  
  #run model function with varied parameter values calling
  source(here::here("script", "07_lhs_model.R"))
  #source("07_lhs_model.R")
  
  #save each run seperately
  saved_carr_dynamics[[nn]] <-  as_data_frame(fitmodel()) %>% dplyr::mutate(run=nn)
  print(paste0("iteration ", nn, " has been done!"))
}

# Save all the runs together
rio::export(saved_carr_dynamics, here::here("results", "saved_carr_dynamicsX.rds"))
