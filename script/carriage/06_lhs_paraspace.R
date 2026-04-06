# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission
# 01/09/2024

#====================================================================

#set seed using a task call for entire session to ensure reproducibility
addTaskCallback(function(...) {set.seed(1988); TRUE})

#set up the parameter space
N_sample_size <-5000

#number of parameters to estimate
N_varied_parameters <- 12

#objects to store sampled probabilities under LHS
uniform_LHS <- randomLHS(N_sample_size, N_varied_parameters)
lhs.sample.mt <- matrix(NA, nrow = N_sample_size, ncol = N_varied_parameters)

#prior uniform distribution of transmission parameters
lhs.sample.mt[,1]  <- qunif(uniform_LHS[,1], min = 1, max = 100) #p1v
lhs.sample.mt[,2]  <- qunif(uniform_LHS[,2], min = 1, max = 100) #p2v
lhs.sample.mt[,3]  <- qunif(uniform_LHS[,3], min = 1, max = 100) #p3v
lhs.sample.mt[,4]  <- qunif(uniform_LHS[,4], min = 1, max = 100) #p4v

lhs.sample.mt[,5]  <- qunif(uniform_LHS[,5], min = 1, max = 100) #p1f
lhs.sample.mt[,6]  <- qunif(uniform_LHS[,6], min = 1, max = 100) #p2f
lhs.sample.mt[,7]  <- qunif(uniform_LHS[,7], min = 1, max = 100) #p3f
lhs.sample.mt[,8]  <- qunif(uniform_LHS[,8], min = 1, max = 100) #p4f

lhs.sample.mt[,9]   <- qunif(uniform_LHS[,9],  min = 1, max = 100) #p1n
lhs.sample.mt[,10]  <- qunif(uniform_LHS[,10], min = 1, max = 100) #p2n
lhs.sample.mt[,11]  <- qunif(uniform_LHS[,11], min = 1, max = 100) #p3n
lhs.sample.mt[,12]  <- qunif(uniform_LHS[,12], min = 1, max = 100) #p4n

# lhs.sample.mt[,1]  <- qunif(uniform_LHS[,1], min = 30, max = 55) #p1v
# lhs.sample.mt[,2]  <- qunif(uniform_LHS[,2], min = 1, max = 2) #p2v
# lhs.sample.mt[,3]  <- qunif(uniform_LHS[,3], min =  0, max = 1) #p3v
# lhs.sample.mt[,4]  <- qunif(uniform_LHS[,4], min =  0, max = 0.5) #p4v

# lhs.sample.mt[,5]  <- qunif(uniform_LHS[,5], min = 0, max = 2) #v1
# lhs.sample.mt[,6]  <- qunif(uniform_LHS[,6], min = 0, max = 2) #v2
# lhs.sample.mt[,7]  <- qunif(uniform_LHS[,7], min = 0, max = 2) #v3
# lhs.sample.mt[,8]  <- qunif(uniform_LHS[,8], min = 0, max = 2) #v4
# 
# lhs.sample.mt[,9]  <- qunif(uniform_LHS[,9], min = 2, max = 3) #f1
# lhs.sample.mt[,10] <- qunif(uniform_LHS[,10], min = 2, max = 3) #f2
# lhs.sample.mt[,11] <- qunif(uniform_LHS[,11], min = 2, max = 3) #f3
# lhs.sample.mt[,12] <- qunif(uniform_LHS[,12], min = 2, max = 3) #f4
# 
# lhs.sample.mt[,13] <- qunif(uniform_LHS[,13], min = 2, max = 3) #n1
# lhs.sample.mt[,14] <- qunif(uniform_LHS[,14], min = 2, max = 3) #n2
# lhs.sample.mt[,15] <- qunif(uniform_LHS[,15], min = 2, max = 3) #n3
# lhs.sample.mt[,16] <- qunif(uniform_LHS[,16], min = 2, max = 3) #n4

#total runs and number of simulations per run
total_run <- 1
len_each_run <- N_sample_size/total_run

rio::export(lhs.sample.mt, file = here("results", "lhs_parsetX.csv"))
