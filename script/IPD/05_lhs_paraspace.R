# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission dynamics and impact of vaccination policies
# 01/09/2024

#====================================================================

#set seed using a task call for entire session to ensure reproducibility
addTaskCallback(function(...) {set.seed(1988); TRUE})

#set up the parameter space
N_sample_size <- 500

#number of parameters to estimate
N_varied_parameters <- 11

#objects to store sampled probabilities under LHS
uniform_LHS <- randomLHS(N_sample_size, N_varied_parameters)
lhs.sample.mt <- matrix(NA, nrow = N_sample_size, ncol = N_varied_parameters)

lhs.sample.mt[,1] <- qunif(uniform_LHS[,1], min = 0.5, max = 0.90) #VE against V carriage - delta1
lhs.sample.mt[,2] <- qunif(uniform_LHS[,2], min = 0.5, max = 0.90) #VE against F carriage - delta2

lhs.sample.mt[,3] <- qunif(uniform_LHS[,3], min = 5, max = 20) #duration of protection of V - omega1
lhs.sample.mt[,4] <- qunif(uniform_LHS[,4], min = 5, max = 20) #duration of protection of F - omega2

lhs.sample.mt[,5] <- qunif(uniform_LHS[,5], min = 0.3, max = 0.7) #generate uncertainty for pre-pcv competition parameter values - compV
lhs.sample.mt[,6] <- qunif(uniform_LHS[,6], min = 0.3, max = 0.7) #generate uncertainty for pre-pcv competition parameter values - compF
lhs.sample.mt[,7] <- qunif(uniform_LHS[,7], min = 0.3, max = 0.7) #generate uncertainty for pre-pcv competition parameter values - compN

lhs.sample.mt[,8] <- qunif(uniform_LHS[,8], min = 1, max = 2) #scale replacement parameter for <1y
lhs.sample.mt[,9] <- qunif(uniform_LHS[,9], min = 3, max = 6) #scale replacement parameter for 1-4y
lhs.sample.mt[,10] <- qunif(uniform_LHS[,10], min = 3, max = 6) #scale replacement parameter for 5-17y
lhs.sample.mt[,11] <- qunif(uniform_LHS[,11], min = 60, max = 100) #scale replacement parameter for 18+y

#total runs and number of simulations per run
total_run <- 1
len_each_run <- N_sample_size/total_run

rio::export(lhs.sample.mt, file = here("results", "IPD", "lhs_parset.csv"))

#turn off the task call to reset seed
removeTaskCallback(1)
