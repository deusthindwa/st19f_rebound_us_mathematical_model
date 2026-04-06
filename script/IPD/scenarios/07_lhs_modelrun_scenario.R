# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission dynamics and impact of vaccination policies
# 01/09/2024

#====================================================================

for (vcase in 8:10) {
  

#load R packages
source(here::here("script", "IPD", "scenarios", "00_main_scenario.R"))

#load data to fit model to
source(here::here("script", "IPD", "scenarios", "01_preprocessData_scenario.R"))

#set up model including births, contacts, population
source(here::here("script", "IPD", "scenarios", "02_setUp_scenario.R"))

#set up fixed parameter values
source(here::here("script", "IPD", "scenarios", "03_modelParKnown_scenario.R"))

#run model function
source(here::here("script", "IPD", "scenarios", "04_transModel_scenario.R"))

#call appropriate latin hypercube samples
#scenario=1, model projections based on estimated values pre-pcv13
#scenario=2, model forward simulation of reduced pcv13 VE against F acquisition
#scenario=3, model forward simulation of reduced pcv13 duration of protection against F carriage
#scenario=4, model forward simulation of reduced susceptibility of F to superinfection
#scenario=5, model fit of reduced pcv13 VE against F acquisition
#scenario=6, model fit of reduced pcv13 duration of protection against F carriage
#scenario=7, model fit of reduced susceptibility of F to superinfection
#scenario=8, model simulation of reduced susceptibility of F to superinfection and reduced Efficacy against F
#scenario=9, model simulation of reduced efficacy against PCV7 serotypes during PCV13 era

scenario = vcase
source(here::here("script", "IPD", "scenarios", "05_lhs_paraspace_scenario.R"))

# Create a list to hold parameter estimates
saved_ipd_modelRun <- list()

#number of model runs
N_run <- 1

for( nn in 1:len_each_run){
  
  #parameter sets to run the dynamic model
  parmAdded <- list(
    delta1 = c(lhs.sample.mt[nn,1]),
    delta2 = c(lhs.sample.mt[nn,2]),
    delta3 = c(lhs.sample.mt[nn,3]),
    delta4 = c(lhs.sample.mt[nn,4]),
    
    omega1 = c(lhs.sample.mt[nn,5]),
    omega2 = c(lhs.sample.mt[nn,6]),
    omega3 = c(lhs.sample.mt[nn,7]),
    omega4 = c(lhs.sample.mt[nn,8]),
    
    compVx  = c(lhs.sample.mt[nn,9]),
    compFx  = c(lhs.sample.mt[nn,10]),
    compNx  = c(lhs.sample.mt[nn,11]),
    
    compVy  = c(lhs.sample.mt[nn,12]),
    compFy  = c(lhs.sample.mt[nn,13]),
    compNy  = c(lhs.sample.mt[nn,14]),
  
    replacex = c(lhs.sample.mt[nn,15], lhs.sample.mt[nn,16], lhs.sample.mt[nn,17], lhs.sample.mt[nn,18]),
    replacey = c(lhs.sample.mt[nn,19], lhs.sample.mt[nn,20], lhs.sample.mt[nn,21], lhs.sample.mt[nn,22])
    )
  
  #add 'parms' without parameters to estimate and add the parm_added with parameters to estimate
  #my_parmset has both fixed and estimated' parameters
  my_parmset <- c(parms, parmAdded)
  
  #run model function with varied parameter values calling
  source(here::here("script", "IPD", "scenarios", "06_lhs_model_scenario.R"))
  
  #save each run separately
  saved_ipd_modelRun[[nn]] <-  as_data_frame(fitmodel()) %>% dplyr::mutate(run=nn)
  print(paste0("iteration ", nn, " has been done!"))
}

#save all the runs together
if (scenario == 1){
  rio::export(saved_ipd_modelRun, here::here("results", "IPD", "saved_ipdRun_projection.rds"))
  
} else if (scenario == 2){
  rio::export(saved_ipd_modelRun, here::here("results", "IPD", "saved_ipdRun_project_efficacy.rds"))
  
} else if (scenario == 3){
  rio::export(saved_ipd_modelRun, here::here("results", "IPD", "saved_ipdRun_project_duration..rds"))
  
} else if (scenario == 4){
  rio::export(saved_ipd_modelRun, here::here("results", "IPD", "saved_ipdRun_project_competition.rds"))
  
} else if (scenario == 5){
  rio::export(saved_ipd_modelRun, here::here("results", "IPD", "saved_ipdRun_modefit_efficacy.rds"))
  
} else if (scenario == 6){
  rio::export(saved_ipd_modelRun, here::here("results", "IPD", "saved_ipdRun_modefit_duration.rds"))
  
} else if (scenario == 7){
  rio::export(saved_ipd_modelRun, here::here("results", "IPD", "saved_ipdRun_modefit_competition.rds"))
  
} else if (scenario == 8){
  rio::export(saved_ipd_modelRun, here::here("results", "IPD", "saved_ipdRun_project_efficacy_reduced.rds"))
  
} else if (scenario == 9){
  rio::export(saved_ipd_modelRun, here::here("results", "IPD", "saved_ipdRun_project_efficacy_constant.rds"))
  
} else if (scenario == 10){
  rio::export(saved_ipd_modelRun, here::here("results", "IPD", "saved_ipdRun_project_pcv7serotypes.rds"))
  
} else{
  print("No save ipdRun")
}

}
