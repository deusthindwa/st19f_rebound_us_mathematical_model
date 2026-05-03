# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission dynamics and impact of vaccination policies
# 01/09/2024

#====================================================================

#load R packages
source(here::here("script", "IPD", "00_main.R"))

#load data to fit model to
source(here::here("script", "IPD", "01_preprocessData.R"))

#set up model including births, contacts, population
source(here::here("script", "IPD", "02_setUp.R"))

#set up fixed parameter values
source(here::here("script", "IPD", "03_modelParKnown.R"))

#run model function
source(here::here("script", "IPD", "04_transModel.R"))

#call latin hypercube samples
source(here::here("script", "IPD", "05_lhs_paraspace.R"))

# Create a list to hold parameter estimates
saved_ipd_dynamics <- list()

#number of model runs
N_run <- 1

for( nn in 1:len_each_run){
  
  #parameter sets to run the dynamic model
  parmAdded <- list(delta1 = c(lhs.sample.mt[nn,1]),
                    delta2 = c(lhs.sample.mt[nn,2]),
                    omega1 = c(lhs.sample.mt[nn,3]),
                    omega2 = c(lhs.sample.mt[nn,4]),
                    compV  = c(lhs.sample.mt[nn,5]),
                    compF  = c(lhs.sample.mt[nn,6]),
                    compN  = c(lhs.sample.mt[nn,7]),
                    replace  = c(lhs.sample.mt[nn,8], lhs.sample.mt[nn,9], lhs.sample.mt[nn,10], lhs.sample.mt[nn,11])
                    )
  
  #add 'parms' without parameters to estimate and add the parm_added with parameters to estimate
  #my_parmset has both fixed and estimated' parameters
  my_parmset <- c(parms, parmAdded)
  
  #run model function with varied parameter values calling
  source(here::here("script", "IPD", "06_lhs_model.R"))
  
  #save each run seperately
  saved_ipd_dynamics[[nn]] <-  as_data_frame(fitmodel()) %>% dplyr::mutate(run=nn)
  print(paste0("iteration ", nn, " has been done!"))
}

#save all the runs together
rio::export(saved_ipd_dynamics, here::here("results", "IPD", "saved_ipd_dynamics.rds"))
