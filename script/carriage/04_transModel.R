# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission
# 01/09/2024

#====================================================================

#source function of the transmission dynamic model
simple_model <- function(t, y, parms, time.step = 'year'){
  
  #read in initial states and their names
  States <- array(y, dim = dim(parms$yinit.matrix))
  dimnames(States) <- dimnames(parms$yinit.matrix)
  
  #unify the time unit of parameter inputs (switch between year/month/week for given parameter values in days)
  if (parms$time.step == 'year') {
    period = 1
    length.step = 365.25 #days
  }
  else if (parms$time.step == 'month') {
    period = 12
    length.step = 30.44 #days
  }
  else { #if (parms$time.step == 'week')
    period = 52.18
    length.step = 7 #days
  }
  
  #pull out the states for the model as vectors
  S <-  States[,'S']
  V <-  States[,'V']
  F <-  States[,'F']
  N <-  States[,'N']
  NV <-  States[,'NV']
  NF <-  States[,'NF']
  VF <-  States[,'VF']
  
  #pull out the number of age groups
  N.ages <- parms$N_ages
  
  #ageing rate (by time step)
  if (parms$time.step == 'year') {
    ageing = 1/(parms$agewidth*period)
    ageout = ageing
    agein = c(1/sum(1/ageing), ageing[-N.ages])
  }
  else if (parms$time.step == 'month') {
    ageing = 1/(parms$agewidth*period)
    ageout = ageing
    agein = c(1/sum(1/ageing), ageing[-N.ages])
  }
  else {#if (parms$time.step == 'week')
    ageing = 1/(parms$agewidth*period)
    ageout = ageing
    agein = c(1/sum(1/ageing), ageing[-N.ages])
  }
  
  #clearance rates of VT, 19F and NVT carriage (r)
  clearV = 1/(parms$durV/length.step)
  clearF = 1/(parms$durF/length.step)
  clearN = 1/(parms$durN/length.step)
  
  #competition parameters
  compV = parms$compV
  compF = parms$compF
  compN = parms$compN
  
  #English population
  pop = df_EW$Population
  
# #============================transmission option 1============================
# 
#   #annual mean number of contacts per person
#   #normalise the contact matrix by dividing by sum of diagonal values
#   #contact_EW = (parms$contact_EW/sum(diag(parms$contact_EW)))*length.step
#   contact_EW = parms$contact_EW*length.step
# 
#   #transmissibility/susceptibility (probability per year)
#   baseV = parms$baseV/(parms$durV/length.step)
#   baseF = parms$baseF/(parms$durF/length.step)
#   baseN = parms$baseN/(parms$durN/length.step)
# 
#   #combines transmissibility or susceptibility per year AND annual mean number of contacts per person
#   #transmissibility=2, for every column, operation | susceptibility=1, for every row, operation
#   #This is the transmission rate per person per year
#   betaV = sweep(contact_EW, 1, baseV, "*")
#   betaF = sweep(contact_EW, 1, baseF, "*")
#   betaN = sweep(contact_EW, 1, baseN, "*")
# 
#   #make transmission rate frequency or density dependent
#   #where (q=0, density) and (q=1, frequency)
#   q = 1
#   betaV = betaV/(sum(yinit.matrix)^(1-q))
#   betaF = betaF/(sum(yinit.matrix)^(1-q))
#   betaN = betaN/(sum(yinit.matrix)^(1-q))
# 
#   #total number of carriers from each ST-group
#   carryV <- (V+VF+NV)/sum(States)*pop
#   carryF <- (F+VF+NF)/sum(States)*pop
#   carryN <- (N+NF+NV)/sum(States)*pop
# 
#   #sum across the columns of contact matrix (beta*I) to compute the force of infection
#   lambdaV <- betaV %*% carryV
#   lambdaF <- betaF %*% carryF
#   lambdaN <- betaN %*% carryN
# 
#   #force of infection plus vectorize force of infection
#   lambdaV <- as.vector(lambdaV)
#   lambdaF <- as.vector(lambdaF)
#   lambdaN <- as.vector(lambdaN)

#============================transmission option 2============================

  #transmission probability per contact
  baseV = parms$baseV
  baseF = parms$baseF
  baseN = parms$baseN

  #transmission probability per contact per unit time
  bV <- baseV / (parms$durV/length.step)
  bF <- baseF / (parms$durF/length.step)
  bN <- baseN / (parms$durN/length.step)

  #let q be transmission type (=0 when density or =1 when frequency-dependent)
  q = 1

  #annual mean number of contacts (from daily contacts)
  contact_EW = parms$contact_EW*length.step

  # Transmission probability per unit time in each age group
  betaV <-  (bV/1000000000)/(sum(yinit.matrix)^(1-q))*contact_EW
  betaF <-  (bF/1000000000)/(sum(yinit.matrix)^(1-q))*contact_EW
  betaN <-  (bN/1000000000)/(sum(yinit.matrix)^(1-q))*contact_EW

  #total number of carriers from each ST-group
  carryV <- (V+VF+NV)/sum(States)*pop
  carryF <- (F+VF+NF)/sum(States)*pop
  carryN <- (N+NF+NV)/sum(States)*pop

  #sum across the columns of contact matrix (beta*I) to compute the force of infection
  lambdaV <- betaV %*% carryV
  lambdaF <- betaF %*% carryF
  lambdaN <- betaN %*% carryN

  #force of infection plus vectorize force of infection
  lambdaV <- as.vector(lambdaV)
  lambdaF <- as.vector(lambdaF)
  lambdaN <- as.vector(lambdaN)
  
  #============================simulate the dynamic model============================
  
  #create a matrix to record the changing variables
  dy <- matrix(NA, nrow = N_ages, ncol = ncol(States))
  colnames(dy) <- colnames(States)
  
  #ordinary differential equations
  dy[,'S'] <- clearV*V + clearN*N + clearF*F -lambdaV*S - lambdaN*S - lambdaF*S - ageout*S + agein*c(sum(States), S[-N.ages])
  
  dy[,'V'] <- lambdaV*S + clearN*NV + clearF*VF - compV*lambdaN*V - compV*lambdaF*V - clearV*V - ageout*V + agein*c(0,V[-N.ages])
  
  dy[,'F'] <- lambdaF*S + clearN*NF + clearV*VF - compF*lambdaN*F - compF*lambdaV*F - clearF*F -  ageout*F + agein*c(0,F[-N.ages])
  
  dy[,'N'] <- lambdaN*S + clearF*NF + clearV*NV - compN*lambdaV*N - compN*lambdaF*N - clearN*N  - ageout*N + agein*c(0,N[-N.ages])
  
  dy[,'NV'] <- compV*lambdaN*V + compN*lambdaV*N - clearV*NV - clearN*NV - ageout*NV + agein*c(0,NV[-N.ages])
  
  dy[,'NF'] <- compN*lambdaF*N + compF*lambdaN*F - clearF*NF - clearN*NF - ageout*NF + agein*c(0,NF[-N.ages])
  
  dy[,'VF'] <- compV*lambdaF*V + compF*lambdaV*F - clearV*VF - clearF*VF - ageout*VF + agein*c(0,VF[-N.ages])
  
  derivs <- as.vector(dy)
  res <- list(derivs)
  
  return(res)
}

#time step is in years
start_time = 1
tmax =  nrow(birth_EW)
my_times <- seq(start_time, tmax, by = 1)


#CHECK PREVALENCE SIMULATIONS BEFORE ANY MODEL FIT

# #using the ODE function to get simulated transmission dynamic results
# results <- ode(y = yinit.vector,
#                t = my_times,
#                func = simple_model,
#                parms = parms)
# plot(results)
# 
# #total population check
# popocheck <- 
#   data.frame(results) %>%
#   dplyr::left_join(data.frame(results) %>% dplyr::mutate(tot = rowSums(.[2:ncol(.)]))) %>%
#   dplyr::select(time, tot)
# 
# popocheck %>%
#   ggplot(aes(x=time, y=tot)) +
#   geom_line(size = 1.5) +
#   theme_bw()
# 
# #total population check
# pop_check <- data.frame(results) %>%
#   left_join(data.frame(results) %>%
#               dplyr::select(everything(), -time) %>%
#               mutate(tot = rowSums(.)))
# 
# #model extraction and prevalence
# res_burned <- results[nrow(results), , drop = FALSE]
# S <- res_burned[,grep('\\b(S)\\b', colnames(res_burned))]
# V <- res_burned[,grep('\\b(V)\\b', colnames(res_burned))]
# F <- res_burned[,grep('\\b(F)\\b', colnames(res_burned))]
# N <- res_burned[,grep('\\b(N)\\b', colnames(res_burned))]
# NV <- res_burned[,grep('\\b(NV)\\b', colnames(res_burned))]
# NF <- res_burned[,grep('\\b(NF)\\b', colnames(res_burned))]
# VF <- res_burned[,grep('\\b(VF)\\b', colnames(res_burned))]
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
#              N.prev.model = Pn)
