# # by Deus Thindwa
# # age-structured mathematical model for pneumococcal transmission
# # 01/09/2024
# 
# #====================================================================
# 
# #set of fixed parameters without parameters to estimate
# parmset <- list(
#   durV = durV,
#   durF = durF,
#   durN = durN,
#   compV = compV,
#   compF = compF,
#   compN = compN,
#   contact_EW = contact_EW,
#   agewidth = agewidth,
#   yinit.matrix = yinit.matrix,
#   N_ages = N_ages,
#   time.step = 'year'
# )
# 
# #time step is in years
# start_time = 1
# tmax =  nrow(birth_EW)
# my_times <- seq(start_time, tmax, by = 1)
# 
# #====================================================================
# 
# fitmodel <- function(parameters, dat) {
# 
#   #transmission probability to estimate
#   baseV <- parameters[(1:parmset$N_ages)]
#   baseF <- parameters[(1:parmset$N_ages)]
#   baseN <- parameters[(1:parmset$N_ages)]
# 
#   compV <- parameters[(parmset$N_ages+1):(2*parmset$N_ages)]
#   compF <- parameters[(2*parmset$N_ages+1):(3*parmset$N_ages)]
#   compN <- parameters[(3*parmset$N_ages+1):(4*parmset$N_ages)]
# 
#   #we need this (double arrow assignment) for the function and model to recognize these parameters
#   baseV <<- baseV
#   baseF <<- baseF
#   baseN <<- baseN
#   # compV <<- compV
#   # compF <<- compF
#   # compN <<- compN
# 
#   #run transmission model with initial conditions & time steps defined above, and parameter values from function call
#   results <- ode(y = yinit.vector, t = my_times,
#                  func = simple_model,
#                  parms = c(parmset,
#                            list(baseV = baseV,
#                                 baseF = baseF,
#                                 baseN = baseN)
#                                 # compV = compV,
#                                 # compF = compF,
#                                 # compN = compN)
#                            )
#                  )
# 
#   res_burned <- results[nrow(results), , drop = FALSE]
#   S <- res_burned[,grep('\\b(S)\\b', colnames(res_burned))]
#   V <- res_burned[,grep('\\b(V)\\b', colnames(res_burned))]
#   F <- res_burned[,grep('\\b(F)\\b', colnames(res_burned))]
#   N <- res_burned[,grep('\\b(N)\\b', colnames(res_burned))]
#   NV <- res_burned[,grep('\\b(NV)\\b', colnames(res_burned))]
#   NF <- res_burned[,grep('\\b(NF)\\b', colnames(res_burned))]
#   VF <- res_burned[,grep('\\b(VF)\\b', colnames(res_burned))]
# 
# #total age-specific carriers after burnin in each obs data available infection states
# #dual carriage are equally shared between states
# S_fit <- S
# V_fit <- V + 0.5*(NV+VF)
# F_fit <- F + 0.5*(NF+VF)
# N_fit <- N + 0.5*(NV+NF)
# 
# #total age-specific infections after burnin
# totInf <- (V_fit + F_fit + N_fit)
# names(totInf) <- c('agegp1 Inf', 'agegp2 Inf', 'agegp3 Inf', 'agegp4 Inf')
# 
# #total age-specific pop after burnin
# totPop <- (V_fit + F_fit + N_fit + S_fit)
# names(totPop) <- c('agegp1 Tot', 'agegp2 Tot', 'agegp3 Tot', 'agegp4 Tot')
# 
# #simulated prevalence
# Ps = S_fit/totPop
# Pv = V_fit/totPop
# Pf = F_fit/totPop
# Pn = N_fit/totPop
# 
# df_EWm <-
#   data_frame(agegp = c("<1y", "1-4y", "5-17y", "18+y"),
#              V.prev.model = Pv,
#              F.prev.model = Pf,
#              N.prev.model = Pn,
#              S.prev.model = Ps)
# 
# df_EWo <- df_EW
# 
# #multinomial likelihood function
# calc_carr_logLik = function(obs_count_carr, model_prev_carr){
#   LL = sum(obs_count_carr * log(model_prev_carr))
#   if(any(model_prev_carr <= 0)) LL = -1000001
#   return(LL)
# }
# 
# #multinomial likelihood for carriage given V, F, N and S states
# LL_prop = 0
# for (age in 1:N_ages){
#   LL_prop = LL_prop +
#     calc_carr_logLik(as.numeric(c(df_EWo$VT.prev[age], df_EWo$F.prev[age], df_EWo$NVT.prev[age], df_EWo$Uncol[age])) * df_EWo$N.total[age],
#                      as.numeric(c(df_EWm$V.prev.model[age], df_EWm$F.prev.model[age], df_EWm$N.prev.model[age], df_EWm$S.prev.model[age]))
#     )
#   return(LL_prop)
# }
# 
# }
# 
# #====================================================================
# 
# #you can minimize the negative log likelihood in order to actually perform the maximum likelihood estimate of the function you are testing.
# #use Optim function for maximum a posteriori (calibrate outputs of the transmission model to observed status)
# #starting values for parameters
# fitLL <- optim(par = c(c(4.512, 1.5, 0.25, 0.25), c(4.512, 1.5, 0.25, 0.25), c(4.512, 1.5, 0.25, 0.25), c(4.512, 1.5, 0.25, 0.25)), #, rep(1, 4), rep(2.5, 4), rep(2.5, 4)
#                fn = fitmodel,
#                dat = sim,
#                control = list(fnscale = -1, maxit = 3000000, trace = TRUE, REPORT = 500))
# 
# #check if the convergence = 0. if not run for several times to make sure the function did not stuck at local minimum
# fitLL
# exp(fitLL$par) # trans prob given contacts
