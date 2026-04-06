# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission dynamics and impact of vaccination policies
# 01/09/2024

#====================================================================

#set up the parameter space
N_sample_size <- 500

#number of parameters to estimate
N_varied_parameters <- 22

#====================================================================

if (scenario == 1){ #model projections based on pcv7 posterior mean estimates

#set seed using a task call for entire session to ensure reproducibility
addTaskCallback(function(...) {set.seed(1988); TRUE})

#objects to store sampled probabilities under LHS
uniform_LHS <- randomLHS(N_sample_size, N_varied_parameters)
lhs.sample.mt <- matrix(NA, nrow = N_sample_size, ncol = N_varied_parameters)

lhs.sample.mt[,1] <- qunif(uniform_LHS[,1], min = 0.808, max = 0.808) #pcv7 VE against V carriage - delta1, fixed
lhs.sample.mt[,2] <- qunif(uniform_LHS[,2], min = 0.676, max = 0.676) #pcv7 VE against F carriage - delta2, fixed

lhs.sample.mt[,3] <- qunif(uniform_LHS[,3], min = 0.808, max = 0.808) #pcv13 VE against V carriage - delta3, fixed 
lhs.sample.mt[,4] <- qunif(uniform_LHS[,4], min = 0.676, max = 0.676) #pcv13 VE against F carriage - delta4, estimate

lhs.sample.mt[,5] <- qunif(uniform_LHS[,5], min = 14.5, max = 14.5) #pcv7 duration of protection of V - omega1, fixed
lhs.sample.mt[,6] <- qunif(uniform_LHS[,6], min = 13.1, max = 13.1) #pcv7 duration of protection of F - omega2, fixed 

lhs.sample.mt[,7] <- qunif(uniform_LHS[,7], min = 14.5, max = 14.5) #pcv13 duration of protection of V - omega3, fixed
lhs.sample.mt[,8] <- qunif(uniform_LHS[,8], min = 13.1, max = 13.1) #pcv13 duration of protection of F - omega4, fixed

lhs.sample.mt[,9] <- qunif(uniform_LHS[,9], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compVx, fixed
lhs.sample.mt[,10] <- qunif(uniform_LHS[,10], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compFx, fixed
lhs.sample.mt[,11] <- qunif(uniform_LHS[,11], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compNx, fixed

lhs.sample.mt[,12] <- qunif(uniform_LHS[,12], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compVy, fixed
lhs.sample.mt[,13] <- qunif(uniform_LHS[,13], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compFy, fixed
lhs.sample.mt[,14] <- qunif(uniform_LHS[,14], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compNy, fixed

lhs.sample.mt[,15] <- qunif(uniform_LHS[,15], min = 1.55, max = 1.55) #scale replacement parameter for <1y
lhs.sample.mt[,16] <- qunif(uniform_LHS[,16], min = 4.38, max = 4.38) #scale replacement parameter for 1-4y
lhs.sample.mt[,17] <- qunif(uniform_LHS[,17], min = 4.53, max = 4.53) #scale replacement parameter for 5-17y
lhs.sample.mt[,18] <- qunif(uniform_LHS[,18], min = 80.3, max = 80.3) #scale replacement parameter for 8+y

lhs.sample.mt[,19] <- qunif(uniform_LHS[,19], min = 0.5, max = 0.5) #scale replacement parameter for <1y
lhs.sample.mt[,20] <- qunif(uniform_LHS[,20], min = 0.1, max = 0.1) #scale replacement parameter for 1-4y
lhs.sample.mt[,21] <- qunif(uniform_LHS[,21], min = 1.1, max = 1.1) #scale replacement parameter for 5-17y
lhs.sample.mt[,22] <- qunif(uniform_LHS[,22], min = 20, max = 20) #scale replacement parameter for 8+y

#total runs and number of simulations per run
total_run <- 1
len_each_run <- N_sample_size/total_run

rio::export(lhs.sample.mt, file = here("results", "IPD", "lhs_parset_projection.csv"))

#turn off the task call to reset seed
removeTaskCallback(1)

#====================================================================

} else if (scenario == 2){ #simulation of reduced efficacy
  
  #set seed using a task call for entire session to ensure reproducibility
  addTaskCallback(function(...) {set.seed(1988); TRUE})
  
  #objects to store sampled probabilities under LHS
  uniform_LHS <- randomLHS(N_sample_size, N_varied_parameters)
  lhs.sample.mt <- matrix(NA, nrow = N_sample_size, ncol = N_varied_parameters)
  
  lhs.sample.mt[,1] <- qunif(uniform_LHS[,1], min = 0.808, max = 0.808) #pcv7 VE against V carriage - delta1, fixed
  lhs.sample.mt[,2] <- qunif(uniform_LHS[,2], min = 0.676, max = 0.676) #pcv7 VE against F carriage - delta2, fixed
  
  lhs.sample.mt[,3] <- qunif(uniform_LHS[,3], min = 0.808, max = 0.808) #pcv13 VE against V carriage - delta3, fixed 
  lhs.sample.mt[,4] <- qunif(uniform_LHS[,4], min = 0.2, max = 0.2) #pcv13 VE against F carriage - delta4, reduced
  
  lhs.sample.mt[,5] <- qunif(uniform_LHS[,5], min = 14.5, max = 14.5) #pcv7 duration of protection of V - omega1, fixed
  lhs.sample.mt[,6] <- qunif(uniform_LHS[,6], min = 13.1, max = 13.1) #pcv7 duration of protection of F - omega2, fixed 
  
  lhs.sample.mt[,7] <- qunif(uniform_LHS[,7], min = 14.5, max = 14.5) #pcv13 duration of protection of V - omega3, fixed
  lhs.sample.mt[,8] <- qunif(uniform_LHS[,8], min = 13.1, max = 13.1) #pcv13 duration of protection of F - omega4, fixed
  
  lhs.sample.mt[,9] <- qunif(uniform_LHS[,9], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compVx, fixed
  lhs.sample.mt[,10] <- qunif(uniform_LHS[,10], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compFx, fixed
  lhs.sample.mt[,11] <- qunif(uniform_LHS[,11], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compNx, fixed
  
  lhs.sample.mt[,12] <- qunif(uniform_LHS[,12], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compVy, fixed
  lhs.sample.mt[,13] <- qunif(uniform_LHS[,13], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compFy, fixed
  lhs.sample.mt[,14] <- qunif(uniform_LHS[,14], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compNy, fixed
  
  lhs.sample.mt[,15] <- qunif(uniform_LHS[,15], min = 1.55, max = 1.55) #scale replacement parameter for <1y
  lhs.sample.mt[,16] <- qunif(uniform_LHS[,16], min = 4.38, max = 4.38) #scale replacement parameter for 1-4y
  lhs.sample.mt[,17] <- qunif(uniform_LHS[,17], min = 4.53, max = 4.53) #scale replacement parameter for 5-17y
  lhs.sample.mt[,18] <- qunif(uniform_LHS[,18], min = 80.3, max = 80.3) #scale replacement parameter for 8+y
  
  lhs.sample.mt[,19] <- qunif(uniform_LHS[,19], min = 0.5, max = 0.5) #scale replacement parameter for <1y
  lhs.sample.mt[,20] <- qunif(uniform_LHS[,20], min = 0.1, max = 0.1) #scale replacement parameter for 1-4y
  lhs.sample.mt[,21] <- qunif(uniform_LHS[,21], min = 1.1, max = 1.1) #scale replacement parameter for 5-17y
  lhs.sample.mt[,22] <- qunif(uniform_LHS[,22], min = 20, max = 20) #scale replacement parameter for 8+y
  
  #total runs and number of simulations per run
  total_run <- 1
  len_each_run <- N_sample_size/total_run
  
  rio::export(lhs.sample.mt, file = here("results", "IPD", "lhs_parset_project_efficacy.csv"))
  
  #turn off the task call to reset seed
  removeTaskCallback(1)
  
  #====================================================================

} else if (scenario == 3){ #simulation of reduced duration of protection
  
  #set seed using a task call for entire session to ensure reproducibility
  addTaskCallback(function(...) {set.seed(1988); TRUE})
  
  #objects to store sampled probabilities under LHS
  uniform_LHS <- randomLHS(N_sample_size, N_varied_parameters)
  lhs.sample.mt <- matrix(NA, nrow = N_sample_size, ncol = N_varied_parameters)
  
  lhs.sample.mt[,1] <- qunif(uniform_LHS[,1], min = 0.808, max = 0.808) #pcv7 VE against V carriage - delta1, fixed
  lhs.sample.mt[,2] <- qunif(uniform_LHS[,2], min = 0.676, max = 0.676) #pcv7 VE against F carriage - delta2, fixed
  
  lhs.sample.mt[,3] <- qunif(uniform_LHS[,3], min = 0.808, max = 0.808) #pcv13 VE against V carriage - delta3, fixed 
  lhs.sample.mt[,4] <- qunif(uniform_LHS[,4], min = 0.676, max = 0.676) #pcv13 VE against F carriage - delta4, fixed
  
  lhs.sample.mt[,5] <- qunif(uniform_LHS[,5], min = 14.5, max = 14.5) #pcv7 duration of protection of V - omega1, fixed
  lhs.sample.mt[,6] <- qunif(uniform_LHS[,6], min = 13.1, max = 13.1) #pcv7 duration of protection of F - omega2, fixed 
  
  lhs.sample.mt[,7] <- qunif(uniform_LHS[,7], min = 14.5, max = 14.5) #pcv13 duration of protection of V - omega3, fixed
  lhs.sample.mt[,8] <- qunif(uniform_LHS[,8], min = 5, max = 5) #pcv13 duration of protection of F - omega4, reduced
  
  lhs.sample.mt[,9] <- qunif(uniform_LHS[,9], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compVx, fixed
  lhs.sample.mt[,10] <- qunif(uniform_LHS[,10], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compFx, fixed
  lhs.sample.mt[,11] <- qunif(uniform_LHS[,11], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compNx, fixed
  
  lhs.sample.mt[,12] <- qunif(uniform_LHS[,12], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compVy, fixed
  lhs.sample.mt[,13] <- qunif(uniform_LHS[,13], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compFy, fixed
  lhs.sample.mt[,14] <- qunif(uniform_LHS[,14], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compNy, fixed
  
  lhs.sample.mt[,15] <- qunif(uniform_LHS[,15], min = 1.55, max = 1.55) #scale replacement parameter for <1y
  lhs.sample.mt[,16] <- qunif(uniform_LHS[,16], min = 4.38, max = 4.38) #scale replacement parameter for 1-4y
  lhs.sample.mt[,17] <- qunif(uniform_LHS[,17], min = 4.53, max = 4.53) #scale replacement parameter for 5-17y
  lhs.sample.mt[,18] <- qunif(uniform_LHS[,18], min = 80.3, max = 80.3) #scale replacement parameter for 8+y
  
  lhs.sample.mt[,19] <- qunif(uniform_LHS[,19], min = 0.5, max = 0.5) #scale replacement parameter for <1y
  lhs.sample.mt[,20] <- qunif(uniform_LHS[,20], min = 0.1, max = 0.1) #scale replacement parameter for 1-4y
  lhs.sample.mt[,21] <- qunif(uniform_LHS[,21], min = 1.1, max = 1.1) #scale replacement parameter for 5-17y
  lhs.sample.mt[,22] <- qunif(uniform_LHS[,22], min = 20, max = 20) #scale replacement parameter for 8+y
  
  #total runs and number of simulations per run
  total_run <- 1
  len_each_run <- N_sample_size/total_run
  
  rio::export(lhs.sample.mt, file = here("results", "IPD", "lhs_parset_project_duration.csv"))
  
  #turn off the task call to reset seed
  removeTaskCallback(1)
  
  #====================================================================
  
} else if (scenario == 4){ #simulation of reduced susceptibility of F 
  
  #set seed using a task call for entire session to ensure reproducibility
  addTaskCallback(function(...) {set.seed(1988); TRUE})
  
  #objects to store sampled probabilities under LHS
  uniform_LHS <- randomLHS(N_sample_size, N_varied_parameters)
  lhs.sample.mt <- matrix(NA, nrow = N_sample_size, ncol = N_varied_parameters)
  
  lhs.sample.mt[,1] <- qunif(uniform_LHS[,1], min = 0.808, max = 0.808) #pcv7 VE against V carriage - delta1, fixed
  lhs.sample.mt[,2] <- qunif(uniform_LHS[,2], min = 0.676, max = 0.676) #pcv7 VE against F carriage - delta2, fixed
  
  lhs.sample.mt[,3] <- qunif(uniform_LHS[,3], min = 0.808, max = 0.808) #pcv13 VE against V carriage - delta3, fixed 
  lhs.sample.mt[,4] <- qunif(uniform_LHS[,4], min = 0.676, max = 0.676) #pcv13 VE against F carriage - delta4, fixed
  
  lhs.sample.mt[,5] <- qunif(uniform_LHS[,5], min = 14.5, max = 14.5) #pcv7 duration of protection of V - omega1, fixed
  lhs.sample.mt[,6] <- qunif(uniform_LHS[,6], min = 13.1, max = 13.1) #pcv7 duration of protection of F - omega2, fixed 
  
  lhs.sample.mt[,7] <- qunif(uniform_LHS[,7], min = 14.5, max = 14.5) #pcv13 duration of protection of V - omega3, fixed
  lhs.sample.mt[,8] <- qunif(uniform_LHS[,8], min = 13.1, max = 13.1) #pcv13 duration of protection of F - omega4, fixed
  
  lhs.sample.mt[,9] <- qunif(uniform_LHS[,9], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compVx, fixed
  lhs.sample.mt[,10] <- qunif(uniform_LHS[,10], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compFx, fixed
  lhs.sample.mt[,11] <- qunif(uniform_LHS[,11], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compNx, fixed
  
  lhs.sample.mt[,12] <- qunif(uniform_LHS[,12], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compVy, fixed
  lhs.sample.mt[,13] <- qunif(uniform_LHS[,13], min = 0.05, max = 0.05) #uncertainty for post-pcv13 competition - compFy, reduced
  lhs.sample.mt[,14] <- qunif(uniform_LHS[,14], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compNy, fixed
  
  lhs.sample.mt[,15] <- qunif(uniform_LHS[,15], min = 1.55, max = 1.55) #scale replacement parameter for <1y
  lhs.sample.mt[,16] <- qunif(uniform_LHS[,16], min = 4.38, max = 4.38) #scale replacement parameter for 1-4y
  lhs.sample.mt[,17] <- qunif(uniform_LHS[,17], min = 4.53, max = 4.53) #scale replacement parameter for 5-17y
  lhs.sample.mt[,18] <- qunif(uniform_LHS[,18], min = 80.3, max = 80.3) #scale replacement parameter for 8+y
  
  lhs.sample.mt[,19] <- qunif(uniform_LHS[,19], min = 0.5, max = 0.5) #scale replacement parameter for <1y
  lhs.sample.mt[,20] <- qunif(uniform_LHS[,20], min = 0.1, max = 0.1) #scale replacement parameter for 1-4y
  lhs.sample.mt[,21] <- qunif(uniform_LHS[,21], min = 1.1, max = 1.1) #scale replacement parameter for 5-17y
  lhs.sample.mt[,22] <- qunif(uniform_LHS[,22], min = 20, max = 20) #scale replacement parameter for 8+y
  
  #total runs and number of simulations per run
  total_run <- 1
  len_each_run <- N_sample_size/total_run
  
  rio::export(lhs.sample.mt, file = here("results", "IPD", "lhs_parset_project_competition.csv"))
  
  #turn off the task call to reset seed
  removeTaskCallback(1)
 
  
  #==================================================================== 
  #====================================================================
  
  
} else if (scenario == 5){ #model fit of reduced efficacy

  #set seed using a task call for entire session to ensure reproducibility
  addTaskCallback(function(...) {set.seed(1988); TRUE})
  
  #objects to store sampled probabilities under LHS
  uniform_LHS <- randomLHS(N_sample_size, N_varied_parameters)
  lhs.sample.mt <- matrix(NA, nrow = N_sample_size, ncol = N_varied_parameters)
  
  lhs.sample.mt[,1] <- qunif(uniform_LHS[,1], min = 0.5, max = 0.9) #pcv7 VE against V carriage - delta1, fixed
  lhs.sample.mt[,2] <- qunif(uniform_LHS[,2], min = 0.5, max = 0.9) #pcv7 VE against F carriage - delta2, fixed
  
  lhs.sample.mt[,3] <- qunif(uniform_LHS[,3], min = 0.5, max = 0.9) #pcv13 VE against V carriage - delta3, fixed 
  lhs.sample.mt[,4] <- qunif(uniform_LHS[,4], min = 0.15, max = 0.25) #pcv13 VE against F carriage - delta4, estimate
  
  lhs.sample.mt[,5] <- qunif(uniform_LHS[,5], min = 5, max = 20) #pcv7 duration of protection of V - omega1, fixed
  lhs.sample.mt[,6] <- qunif(uniform_LHS[,6], min = 5, max = 20) #pcv7 duration of protection of F - omega2, fixed 
  
  lhs.sample.mt[,7] <- qunif(uniform_LHS[,7], min = 5, max = 20) #pcv13 duration of protection of V - omega3, fixed
  lhs.sample.mt[,8] <- qunif(uniform_LHS[,8], min = 5, max = 20) #pcv13 duration of protection of F - omega4, fixed
  
  lhs.sample.mt[,9] <- qunif(uniform_LHS[,9], min = 0.35, max = 0.45) #uncertainty for pre-pcv13 competition - compVx, fixed
  lhs.sample.mt[,10] <- qunif(uniform_LHS[,10], min = 0.35, max = 0.45) #uncertainty for pre-pcv13 competition - compFx, fixed
  lhs.sample.mt[,11] <- qunif(uniform_LHS[,11], min = 0.35, max = 0.45) #uncertainty for pre-pcv13 competition - compNx, fixed
  
  lhs.sample.mt[,12] <- qunif(uniform_LHS[,12], min = 0.35, max = 0.45) #uncertainty for post-pcv13 competition - compVy, fixed
  lhs.sample.mt[,13] <- qunif(uniform_LHS[,13], min = 0.35, max = 0.45) #uncertainty for post-pcv13 competition - compFy, fixed
  lhs.sample.mt[,14] <- qunif(uniform_LHS[,14], min = 0.35, max = 0.45) #uncertainty for post-pcv13 competition - compNy, fixed
  
  lhs.sample.mt[,15] <- qunif(uniform_LHS[,15], min = 1.55, max = 1.55) #scale replacement parameter for <1y
  lhs.sample.mt[,16] <- qunif(uniform_LHS[,16], min = 4.38, max = 4.38) #scale replacement parameter for 1-4y
  lhs.sample.mt[,17] <- qunif(uniform_LHS[,17], min = 4.53, max = 4.53) #scale replacement parameter for 5-17y
  lhs.sample.mt[,18] <- qunif(uniform_LHS[,18], min = 80.3, max = 80.3) #scale replacement parameter for 8+y
  
  lhs.sample.mt[,19] <- qunif(uniform_LHS[,19], min = 0.5, max = 0.5) #scale replacement parameter for <1y
  lhs.sample.mt[,20] <- qunif(uniform_LHS[,20], min = 0.1, max = 0.1) #scale replacement parameter for 1-4y
  lhs.sample.mt[,21] <- qunif(uniform_LHS[,21], min = 1.1, max = 1.1) #scale replacement parameter for 5-17y
  lhs.sample.mt[,22] <- qunif(uniform_LHS[,22], min = 20, max = 20) #scale replacement parameter for 8+y
  
  #total runs and number of simulations per run
  total_run <- 1
  len_each_run <- N_sample_size/total_run
  
  rio::export(lhs.sample.mt, file = here("results", "IPD", "lhs_parset_modefit_efficacy.csv"))
  
  #turn off the task call to reset seed
  removeTaskCallback(1)
  
  #==================================================================== 
  
} else if (scenario == 6){ #model fit of reduced duration of protection
  
  #set seed using a task call for entire session to ensure reproducibility
  addTaskCallback(function(...) {set.seed(1988); TRUE})
  
  #objects to store sampled probabilities under LHS
  uniform_LHS <- randomLHS(N_sample_size, N_varied_parameters)
  lhs.sample.mt <- matrix(NA, nrow = N_sample_size, ncol = N_varied_parameters)
  
  lhs.sample.mt[,1] <- qunif(uniform_LHS[,1], min = 0.5, max = 0.9) #pcv7 VE against V carriage - delta1, fixed
  lhs.sample.mt[,2] <- qunif(uniform_LHS[,2], min = 0.5, max = 0.9) #pcv7 VE against F carriage - delta2, fixed
  
  lhs.sample.mt[,3] <- qunif(uniform_LHS[,3], min = 0.5, max = 0.9) #pcv13 VE against V carriage - delta3, fixed 
  lhs.sample.mt[,4] <- qunif(uniform_LHS[,4], min = 0.5, max = 0.9) #pcv13 VE against F carriage - delta4, estimate
  
  lhs.sample.mt[,5] <- qunif(uniform_LHS[,5], min = 5, max = 20) #pcv7 duration of protection of V - omega1, fixed
  lhs.sample.mt[,6] <- qunif(uniform_LHS[,6], min = 5, max = 20) #pcv7 duration of protection of F - omega2, fixed 
  
  lhs.sample.mt[,7] <- qunif(uniform_LHS[,7], min = 5, max = 20) #pcv13 duration of protection of V - omega3, fixed
  lhs.sample.mt[,8] <- qunif(uniform_LHS[,8], min = 2, max = 8) #pcv13 duration of protection of F - omega4, fixed
  
  lhs.sample.mt[,9] <- qunif(uniform_LHS[,9], min = 0.35, max = 0.45) #uncertainty for pre-pcv13 competition - compVx, fixed
  lhs.sample.mt[,10] <- qunif(uniform_LHS[,10], min = 0.35, max = 0.45) #uncertainty for pre-pcv13 competition - compFx, fixed
  lhs.sample.mt[,11] <- qunif(uniform_LHS[,11], min = 0.35, max = 0.45) #uncertainty for pre-pcv13 competition - compNx, fixed
  
  lhs.sample.mt[,12] <- qunif(uniform_LHS[,12], min = 0.35, max = 0.45) #uncertainty for post-pcv13 competition - compVy, fixed
  lhs.sample.mt[,13] <- qunif(uniform_LHS[,13], min = 0.35, max = 0.45) #uncertainty for post-pcv13 competition - compFy, fixed
  lhs.sample.mt[,14] <- qunif(uniform_LHS[,14], min = 0.35, max = 0.45) #uncertainty for post-pcv13 competition - compNy, fixed
  
  lhs.sample.mt[,15] <- qunif(uniform_LHS[,15], min = 1.55, max = 1.55) #scale replacement parameter for <1y
  lhs.sample.mt[,16] <- qunif(uniform_LHS[,16], min = 4.38, max = 4.38) #scale replacement parameter for 1-4y
  lhs.sample.mt[,17] <- qunif(uniform_LHS[,17], min = 4.53, max = 4.53) #scale replacement parameter for 5-17y
  lhs.sample.mt[,18] <- qunif(uniform_LHS[,18], min = 80.3, max = 80.3) #scale replacement parameter for 8+y
  
  lhs.sample.mt[,19] <- qunif(uniform_LHS[,19], min = 0.5, max = 0.5) #scale replacement parameter for <1y
  lhs.sample.mt[,20] <- qunif(uniform_LHS[,20], min = 0.1, max = 0.1) #scale replacement parameter for 1-4y
  lhs.sample.mt[,21] <- qunif(uniform_LHS[,21], min = 1.1, max = 1.1) #scale replacement parameter for 5-17y
  lhs.sample.mt[,22] <- qunif(uniform_LHS[,22], min = 20, max = 20) #scale replacement parameter for 8+y
  
  #total runs and number of simulations per run
  total_run <- 1
  len_each_run <- N_sample_size/total_run
  
  rio::export(lhs.sample.mt, file = here("results", "IPD", "lhs_parset_modefit_duration.csv"))
  
  #turn off the task call to reset seed
  removeTaskCallback(1)
  
  #==================================================================== 
  
} else if (scenario == 7) { #model fit of reduced susceptibility of F 
  
  #set seed using a task call for entire session to ensure reproducibility
  addTaskCallback(function(...) {set.seed(1988); TRUE})
  
  #objects to store sampled probabilities under LHS
  uniform_LHS <- randomLHS(N_sample_size, N_varied_parameters)
  lhs.sample.mt <- matrix(NA, nrow = N_sample_size, ncol = N_varied_parameters)
  
  lhs.sample.mt[,1] <- qunif(uniform_LHS[,1], min = 0.5, max = 0.9) #pcv7 VE against V carriage - delta1, fixed
  lhs.sample.mt[,2] <- qunif(uniform_LHS[,2], min = 0.5, max = 0.9) #pcv7 VE against F carriage - delta2, fixed
  
  lhs.sample.mt[,3] <- qunif(uniform_LHS[,3], min = 0.5, max = 0.9) #pcv13 VE against V carriage - delta3, fixed 
  lhs.sample.mt[,4] <- qunif(uniform_LHS[,4], min = 0.5, max = 0.9) #pcv13 VE against F carriage - delta4, estimate
  
  lhs.sample.mt[,5] <- qunif(uniform_LHS[,5], min = 5, max = 20) #pcv7 duration of protection of V - omega1, fixed
  lhs.sample.mt[,6] <- qunif(uniform_LHS[,6], min = 5, max = 20) #pcv7 duration of protection of F - omega2, fixed 
  
  lhs.sample.mt[,7] <- qunif(uniform_LHS[,7], min = 5, max = 20) #pcv13 duration of protection of V - omega3, fixed
  lhs.sample.mt[,8] <- qunif(uniform_LHS[,8], min = 5, max = 20) #pcv13 duration of protection of F - omega4, fixed
  
  lhs.sample.mt[,9] <- qunif(uniform_LHS[,9], min = 0.35, max = 0.45) #uncertainty for pre-pcv13 competition - compVx, fixed
  lhs.sample.mt[,10] <- qunif(uniform_LHS[,10], min = 0.35, max = 0.45) #uncertainty for pre-pcv13 competition - compFx, fixed
  lhs.sample.mt[,11] <- qunif(uniform_LHS[,11], min = 0.35, max = 0.45) #uncertainty for pre-pcv13 competition - compNx, fixed
  
  lhs.sample.mt[,12] <- qunif(uniform_LHS[,12], min = 0.35, max = 0.45) #uncertainty for post-pcv13 competition - compVy, fixed
  lhs.sample.mt[,13] <- qunif(uniform_LHS[,13], min = 0.01, max = 0.1) #uncertainty for post-pcv13 competition - compFy, fixed
  lhs.sample.mt[,14] <- qunif(uniform_LHS[,14], min = 0.35, max = 0.45) #uncertainty for post-pcv13 competition - compNy, fixed
  
  lhs.sample.mt[,15] <- qunif(uniform_LHS[,15], min = 1.55, max = 1.55) #scale replacement parameter for <1y
  lhs.sample.mt[,16] <- qunif(uniform_LHS[,16], min = 4.38, max = 4.38) #scale replacement parameter for 1-4y
  lhs.sample.mt[,17] <- qunif(uniform_LHS[,17], min = 4.53, max = 4.53) #scale replacement parameter for 5-17y
  lhs.sample.mt[,18] <- qunif(uniform_LHS[,18], min = 80.3, max = 80.3) #scale replacement parameter for 8+y
  
  lhs.sample.mt[,19] <- qunif(uniform_LHS[,19], min = 0.5, max = 0.5) #scale replacement parameter for <1y
  lhs.sample.mt[,20] <- qunif(uniform_LHS[,20], min = 0.1, max = 0.1) #scale replacement parameter for 1-4y
  lhs.sample.mt[,21] <- qunif(uniform_LHS[,21], min = 1.1, max = 1.1) #scale replacement parameter for 5-17y
  lhs.sample.mt[,22] <- qunif(uniform_LHS[,22], min = 20, max = 20) #scale replacement parameter for 8+y
  
  #total runs and number of simulations per run
  total_run <- 1
  len_each_run <- N_sample_size/total_run
  
  rio::export(lhs.sample.mt, file = here("results", "IPD", "lhs_parset_modefit_competition.csv"))
  
  #turn off the task call to reset seed
  removeTaskCallback(1)
  
  
  #==================================================================== 
  #==================================================================== 
  
  
} else if (scenario == 8){ #simulation of reduced efficacy for both PCV13 against both V and F
  
  #set seed using a task call for entire session to ensure reproducibility
  addTaskCallback(function(...) {set.seed(1988); TRUE})
  
  #objects to store sampled probabilities under LHS
  uniform_LHS <- randomLHS(N_sample_size, N_varied_parameters)
  lhs.sample.mt <- matrix(NA, nrow = N_sample_size, ncol = N_varied_parameters)
  
  lhs.sample.mt[,1] <- qunif(uniform_LHS[,1], min = 0.808, max = 0.808) #pcv7 VE against V carriage - delta1, fixed
  lhs.sample.mt[,2] <- qunif(uniform_LHS[,2], min = 0.808, max = 0.808) #pcv7 VE against F carriage - delta2, fixed
  
  lhs.sample.mt[,3] <- qunif(uniform_LHS[,3], min = 0.2, max = 0.2) #pcv13 VE against V carriage - delta3, fixed 
  lhs.sample.mt[,4] <- qunif(uniform_LHS[,4], min = 0.2, max = 0.2) #pcv13 VE against F carriage - delta4, fixed
  
  lhs.sample.mt[,5] <- qunif(uniform_LHS[,5], min = 14.5, max = 14.5) #pcv7 duration of protection of V - omega1, fixed
  lhs.sample.mt[,6] <- qunif(uniform_LHS[,6], min = 13.1, max = 13.1) #pcv7 duration of protection of F - omega2, fixed 
  
  lhs.sample.mt[,7] <- qunif(uniform_LHS[,7], min = 14.5, max = 14.5) #pcv13 duration of protection of V - omega3, fixed
  lhs.sample.mt[,8] <- qunif(uniform_LHS[,8], min = 13.1, max = 13.1) #pcv13 duration of protection of F - omega4, fixed
  
  lhs.sample.mt[,9] <- qunif(uniform_LHS[,9], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compVx, fixed
  lhs.sample.mt[,10] <- qunif(uniform_LHS[,10], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compFx, fixed
  lhs.sample.mt[,11] <- qunif(uniform_LHS[,11], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compNx, fixed
  
  lhs.sample.mt[,12] <- qunif(uniform_LHS[,12], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compVy, fixed
  lhs.sample.mt[,13] <- qunif(uniform_LHS[,13], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compFy, reduced
  lhs.sample.mt[,14] <- qunif(uniform_LHS[,14], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compNy, fixed
  
  lhs.sample.mt[,15] <- qunif(uniform_LHS[,15], min = 1.55, max = 1.55) #scale replacement parameter for <1y
  lhs.sample.mt[,16] <- qunif(uniform_LHS[,16], min = 4.38, max = 4.38) #scale replacement parameter for 1-4y
  lhs.sample.mt[,17] <- qunif(uniform_LHS[,17], min = 4.53, max = 4.53) #scale replacement parameter for 5-17y
  lhs.sample.mt[,18] <- qunif(uniform_LHS[,18], min = 80.3, max = 80.3) #scale replacement parameter for 8+y
  
  lhs.sample.mt[,19] <- qunif(uniform_LHS[,19], min = 0.5, max = 0.5) #scale replacement parameter for <1y
  lhs.sample.mt[,20] <- qunif(uniform_LHS[,20], min = 0.1, max = 0.1) #scale replacement parameter for 1-4y
  lhs.sample.mt[,21] <- qunif(uniform_LHS[,21], min = 1.1, max = 1.1) #scale replacement parameter for 5-17y
  lhs.sample.mt[,22] <- qunif(uniform_LHS[,22], min = 20, max = 20) #scale replacement parameter for 8+y
  
  #total runs and number of simulations per run
  total_run <- 1
  len_each_run <- N_sample_size/total_run
  
  rio::export(lhs.sample.mt, file = here("results", "IPD", "lhs_parset_project_efficacy_reduced.csv"))
  
  #turn off the task call to reset seed
  removeTaskCallback(1)
  
  #==================================================================== 

} else if (scenario == 9){ #simulation of constant efficacy for PCV7 and PCV13
  
  #set seed using a task call for entire session to ensure reproducibility
  addTaskCallback(function(...) {set.seed(1988); TRUE})
  
  #objects to store sampled probabilities under LHS
  uniform_LHS <- randomLHS(N_sample_size, N_varied_parameters)
  lhs.sample.mt <- matrix(NA, nrow = N_sample_size, ncol = N_varied_parameters)
  
  lhs.sample.mt[,1] <- qunif(uniform_LHS[,1], min = 0.808, max = 0.808) #pcv7 VE against V carriage - delta1, fixed
  lhs.sample.mt[,2] <- qunif(uniform_LHS[,2], min = 0.808, max = 0.808) #pcv7 VE against F carriage - delta2, fixed
  
  lhs.sample.mt[,3] <- qunif(uniform_LHS[,3], min = 0.808, max = 0.808) #pcv13 VE against V carriage - delta3, fixed 
  lhs.sample.mt[,4] <- qunif(uniform_LHS[,4], min = 0.808, max = 0.808) #pcv13 VE against F carriage - delta4, reduced
  
  lhs.sample.mt[,5] <- qunif(uniform_LHS[,5], min = 14.5, max = 14.5) #pcv7 duration of protection of V - omega1, fixed
  lhs.sample.mt[,6] <- qunif(uniform_LHS[,6], min = 13.1, max = 13.1) #pcv7 duration of protection of F - omega2, fixed 
  
  lhs.sample.mt[,7] <- qunif(uniform_LHS[,7], min = 14.5, max = 14.5) #pcv13 duration of protection of V - omega3, fixed
  lhs.sample.mt[,8] <- qunif(uniform_LHS[,8], min = 13.1, max = 13.1) #pcv13 duration of protection of F - omega4, fixed
  
  lhs.sample.mt[,9] <- qunif(uniform_LHS[,9], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compVx, fixed
  lhs.sample.mt[,10] <- qunif(uniform_LHS[,10], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compFx, fixed
  lhs.sample.mt[,11] <- qunif(uniform_LHS[,11], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compNx, fixed
  
  lhs.sample.mt[,12] <- qunif(uniform_LHS[,12], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compVy, fixed
  lhs.sample.mt[,13] <- qunif(uniform_LHS[,13], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compFy, fixed
  lhs.sample.mt[,14] <- qunif(uniform_LHS[,14], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compNy, fixed
  
  lhs.sample.mt[,15] <- qunif(uniform_LHS[,15], min = 1.55, max = 1.55) #scale replacement parameter for <1y
  lhs.sample.mt[,16] <- qunif(uniform_LHS[,16], min = 4.38, max = 4.38) #scale replacement parameter for 1-4y
  lhs.sample.mt[,17] <- qunif(uniform_LHS[,17], min = 4.53, max = 4.53) #scale replacement parameter for 5-17y
  lhs.sample.mt[,18] <- qunif(uniform_LHS[,18], min = 80.3, max = 80.3) #scale replacement parameter for 8+y
  
  lhs.sample.mt[,19] <- qunif(uniform_LHS[,19], min = 0.5, max = 0.5) #scale replacement parameter for <1y
  lhs.sample.mt[,20] <- qunif(uniform_LHS[,20], min = 0.1, max = 0.1) #scale replacement parameter for 1-4y
  lhs.sample.mt[,21] <- qunif(uniform_LHS[,21], min = 1.1, max = 1.1) #scale replacement parameter for 5-17y
  lhs.sample.mt[,22] <- qunif(uniform_LHS[,22], min = 20, max = 20) #scale replacement parameter for 8+y
  
  #total runs and number of simulations per run
  total_run <- 1
  len_each_run <- N_sample_size/total_run
  
  rio::export(lhs.sample.mt, file = here("results", "IPD", "lhs_parset_project_efficacy_constant.csv"))
  
  #turn off the task call to reset seed
  removeTaskCallback(1)
  
} else if (scenario == 10){ #simulation of reduced efficacy for both PCV13 against both V and F but projection for PCV7 ST only
  
  #set seed using a task call for entire session to ensure reproducibility
  addTaskCallback(function(...) {set.seed(1988); TRUE})
  
  #objects to store sampled probabilities under LHS
  uniform_LHS <- randomLHS(N_sample_size, N_varied_parameters)
  lhs.sample.mt <- matrix(NA, nrow = N_sample_size, ncol = N_varied_parameters)
  
  lhs.sample.mt[,1] <- qunif(uniform_LHS[,1], min = 0.808, max = 0.808) #pcv7 VE against V carriage - delta1, fixed
  lhs.sample.mt[,2] <- qunif(uniform_LHS[,2], min = 0.808, max = 0.808) #pcv7 VE against F carriage - delta2, fixed
  
  lhs.sample.mt[,3] <- qunif(uniform_LHS[,3], min = 0.2, max = 0.2) #pcv13 VE against V carriage - delta3, fixed 
  lhs.sample.mt[,4] <- qunif(uniform_LHS[,4], min = 0.2, max = 0.2) #pcv13 VE against F carriage - delta4, fixed
  
  lhs.sample.mt[,5] <- qunif(uniform_LHS[,5], min = 14.5, max = 14.5) #pcv7 duration of protection of V - omega1, fixed
  lhs.sample.mt[,6] <- qunif(uniform_LHS[,6], min = 13.1, max = 13.1) #pcv7 duration of protection of F - omega2, fixed 
  
  lhs.sample.mt[,7] <- qunif(uniform_LHS[,7], min = 14.5, max = 14.5) #pcv13 duration of protection of V - omega3, fixed
  lhs.sample.mt[,8] <- qunif(uniform_LHS[,8], min = 13.1, max = 13.1) #pcv13 duration of protection of F - omega4, fixed
  
  lhs.sample.mt[,9] <- qunif(uniform_LHS[,9], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compVx, fixed
  lhs.sample.mt[,10] <- qunif(uniform_LHS[,10], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compFx, fixed
  lhs.sample.mt[,11] <- qunif(uniform_LHS[,11], min = 0.4, max = 0.4) #uncertainty for pre-pcv13 competition - compNx, fixed
  
  lhs.sample.mt[,12] <- qunif(uniform_LHS[,12], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compVy, fixed
  lhs.sample.mt[,13] <- qunif(uniform_LHS[,13], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compFy, reduced
  lhs.sample.mt[,14] <- qunif(uniform_LHS[,14], min = 0.4, max = 0.4) #uncertainty for post-pcv13 competition - compNy, fixed
  
  lhs.sample.mt[,15] <- qunif(uniform_LHS[,15], min = 1.55, max = 1.55) #scale replacement parameter for <1y
  lhs.sample.mt[,16] <- qunif(uniform_LHS[,16], min = 4.38, max = 4.38) #scale replacement parameter for 1-4y
  lhs.sample.mt[,17] <- qunif(uniform_LHS[,17], min = 4.53, max = 4.53) #scale replacement parameter for 5-17y
  lhs.sample.mt[,18] <- qunif(uniform_LHS[,18], min = 80.3, max = 80.3) #scale replacement parameter for 8+y
  
  lhs.sample.mt[,19] <- qunif(uniform_LHS[,19], min = 0.5, max = 0.5) #scale replacement parameter for <1y
  lhs.sample.mt[,20] <- qunif(uniform_LHS[,20], min = 0.1, max = 0.1) #scale replacement parameter for 1-4y
  lhs.sample.mt[,21] <- qunif(uniform_LHS[,21], min = 1.1, max = 1.1) #scale replacement parameter for 5-17y
  lhs.sample.mt[,22] <- qunif(uniform_LHS[,22], min = 20, max = 20) #scale replacement parameter for 8+y
  
  #total runs and number of simulations per run
  total_run <- 1
  len_each_run <- N_sample_size/total_run
  
  rio::export(lhs.sample.mt, file = here("results", "IPD", "lhs_parset_project_pcv7serotypes.csv"))
  
  #turn off the task call to reset seed
  removeTaskCallback(1)
  
} else {
  
  print("No scenario selected")
  
}
